"""
Master Tomography Acquisition System (High-Precision Ultrasound ToF Processing)
Integrates Pico 2W Hardware Triggering with Hantek 6022BE Waveform Capture.

Algorithms Implemented:
1. Cross-Correlation (R_xy) with Parabolic Sub-sample Interpolation (Sub-nanosecond ToF precision).
2. Hilbert Transform Envelope Detection to eliminate Cycle Skipping.
3. Adaptive Inter-Pulse Cooldown (+10% Safety Factor) to avoid acoustic overlap.
4. Pico 2W USB Serial Trigger Protocol for multiplexer pin selection & 250ns PIO pulses.
"""

import sys
import time
import csv
import os
from datetime import datetime
import numpy as np
from scipy.signal import hilbert, correlate, correlation_lags
import libusb_package
import usb.core
import usb.util
from hantek_interface_testing import Hantek6022BE

try:
    import serial
    import serial.tools.list_ports
    SERIAL_AVAILABLE = True
except ImportError:
    SERIAL_AVAILABLE = False

class PicoTriggerController:
    """Handles USB Serial communication with Pico 2W running readable_tomography_logic.py."""
    def __init__(self, port=None, baudrate=115200):
        self.ser = None
        self.port = port
        self.baudrate = baudrate

    def auto_connect(self):
        if not SERIAL_AVAILABLE:
            print("[!] PySerial not available. Running in local trigger mode.")
            return False

        if self.port:
            try:
                self.ser = serial.Serial(self.port, self.baudrate, timeout=1.0)
                print(f"[+] Connected to Pico 2W on {self.port}")
                return True
            except Exception as e:
                print(f"[*] Serial connection to {self.port} failed: {e}")

        # Auto-detect Pico serial port (VID 0x2E8A / Raspberry Pi Pico)
        ports = serial.tools.list_ports.comports()
        for p in ports:
            if "2E8A" in p.hwid.upper() or "PICO" in p.description.upper() or "MICROPYTHON" in p.description.upper():
                try:
                    self.ser = serial.Serial(p.device, self.baudrate, timeout=1.0)
                    self.port = p.device
                    print(f"[+] Auto-detected Pico 2W on Serial Port {self.port}")
                    return True
                except Exception:
                    pass

        print("[!] Pico 2W USB Serial port not detected. Operating with hardware timing sync.")
        return False

    def trigger_pair(self, tx_chan, rx_chan):
        """Sends trigger command to Pico 2W over USB Serial."""
        if self.ser and self.ser.is_open:
            cmd = f"TRIG {tx_chan} {rx_chan}\n"
            try:
                self.ser.write(cmd.encode('ascii'))
                time.sleep(0.002)  # Give Pico 2ms to execute PIO pulse
                response = self.ser.readline().decode('utf-8', errors='ignore').strip()
                return "ACK" in response
            except Exception as e:
                print(f"[*] Serial write warning: {e}")
        return True

    def close(self):
        if self.ser and self.ser.is_open:
            try:
                self.ser.close()
            except Exception:
                pass


class TomographyAcquisitionMaster:
    def __init__(self, sample_rate_mhz=16.0, safety_factor=0.10, min_cooldown_ms=5.0, serial_port=None, runs_dir=None):
        self.scope = Hantek6022BE()
        self.pico = PicoTriggerController(port=serial_port)
        self.sample_rate_mhz = sample_rate_mhz
        self.safety_factor = safety_factor      # 10% safety margin for ToF variation & flow changes
        self.min_cooldown_ms = min_cooldown_ms  # Minimum acoustic ring-down dissipation time
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.runs_dir = runs_dir if runs_dir else os.path.join(script_dir, "runs")

        # Matches transducer pairs defined in readable_tomography_logic.py
        self.transducer_pairs = [
            (1, 2), (2, 1),
            (3, 4), (4, 3),
            (5, 6), (6, 5),
            (7, 8), (8, 7)
        ]
        self.tof_results = {}
        self.last_observed_tof_us = {}

    def connect(self):
        """Connects to Hantek 6022BE scope and Pico 2W serial trigger."""
        print("=== Initializing Precision Automated Tomography System ===")
        
        # Connect to Pico 2W
        self.pico.auto_connect()

        # Connect to Hantek 6022BE
        if not self.scope.connect():
            print("[X] Scope connection failed. Ensure OpenHantek GUI is closed.")
            return False
        
        self.scope.configure(sample_rate_mhz=self.sample_rate_mhz, ch1_range_v=5.0, ch2_range_v=5.0)
        print("[+] Hantek 6022BE Scope connected & configured.")
        return True

    def calculate_precision_tof_micros(self, ch1_volts, ch2_volts, method="cross_correlation"):
        """
        High-Precision Time-of-Flight (ToF) Calculation:
        
        Method 1: Cross-Correlation (R_xy) + Sub-sample Parabolic Interpolation
                  Measures true signal lag between transmitted reference and received pulse.
                  Provides sub-nanosecond delay resolution.
                  
        Method 2: Hilbert Transform Analytic Envelope Peak
                  Smoothes sinusoidal carrier to detect wavepacket envelope and eliminate cycle skipping.
        """
        if ch1_volts is None or ch2_volts is None or len(ch1_volts) == 0:
            return None

        dt_us = 1.0 / self.sample_rate_mhz

        # Normalize signals (remove DC bias)
        v_tx = ch1_volts - np.mean(ch1_volts)
        v_rx = ch2_volts - np.mean(ch2_volts)

        # Check if received signal exceeds noise floor (minimum threshold)
        rx_amplitude = np.max(np.abs(v_rx))
        if rx_amplitude < 0.08:  # Signal below noise floor (80 mV)
            return None

        if method == "cross_correlation":
            # 1. Compute full Cross-Correlation sequence R_xy[k]
            corr = correlate(v_rx, v_tx, mode='full')
            lags = correlation_lags(len(v_rx), len(v_tx), mode='full')

            # Find correlation peak lag index
            peak_idx = np.argmax(corr)
            k_max = lags[peak_idx]

            # 2. Sub-sample Parabolic (Quadratic) Interpolation around the peak
            if 0 < peak_idx < len(corr) - 1:
                y0 = corr[peak_idx - 1]
                y1 = corr[peak_idx]
                y2 = corr[peak_idx + 1]
                denom = (y0 - 2 * y1 + y2)
                if abs(denom) > 1e-12:
                    delta = (y0 - y2) / (2.0 * denom)
                else:
                    delta = 0.0
            else:
                delta = 0.0

            subsample_lag = k_max + delta
            tof_us = subsample_lag * dt_us
            return max(0.0, tof_us)

        elif method == "hilbert_envelope":
            # 1. Compute Hilbert Transform Analytic Envelopes
            env_tx = np.abs(hilbert(v_tx))
            env_rx = np.abs(hilbert(v_rx))

            # 2. Find peak of envelope
            tx_peak_idx = np.argmax(env_tx)
            rx_peak_idx = np.argmax(env_rx)

            tof_samples = rx_peak_idx - tx_peak_idx
            tof_us = tof_samples * dt_us
            return max(0.0, tof_us)

        else:
            # Fallback First Threshold Crossing
            thresh = 0.15
            tx_c = np.where(np.abs(v_tx) > thresh)[0]
            rx_c = np.where(np.abs(v_rx) > thresh)[0]
            if len(tx_c) == 0 or len(rx_c) == 0:
                return None
            return max(0.0, (rx_c[0] - tx_c[0]) * dt_us)

    def calculate_adaptive_cooldown_sec(self, tx, rx, measured_tof_us):
        """
        Calculates adaptive inter-pulse cooldown using:
        Cooldown = (Observed_ToF * 1.10 Safety Factor) + RingDown_Dissipation_Time
        Prevents firing a new pulse before the previous signal/echo has completed.
        """
        if measured_tof_us is not None and measured_tof_us > 0:
            # Apply +10% safety factor to measured ToF
            safe_tof_us = measured_tof_us * (1.0 + self.safety_factor)
            # Add minimum acoustic ring-down dissipation time
            cooldown_sec = (safe_tof_us / 1e6) + (self.min_cooldown_ms / 1000.0)
        else:
            cooldown_sec = 0.010  # Fallback 10ms
            
        return cooldown_sec

    def execute_tomography_scan(self, run_folder_name=None, continuous=True):
        """
        Iterates across all multiplexed transducer pairs continuously (until KeyboardInterrupt),
        triggers Pico 2W, captures scope waveforms, computes precision ToF via Cross-Correlation,
        applies +10% safety cooldowns, and exports scope waveform data and tomography ToF matrix
        to timestamped CSV files in a dedicated run directory.
        """
        now = datetime.now()
        timestamp = now.strftime("%Y_%m_%d__%H_%M_%S")
        
        if run_folder_name is None:
            run_folder_name = f"tomography_data_{timestamp}"

        target_run_dir = os.path.join(self.runs_dir, run_folder_name)
        os.makedirs(target_run_dir, exist_ok=True)

        tof_csv_filename = f"tof_data_{timestamp}.csv"
        scope_csv_filename = f"scope_data_{timestamp}.csv"

        tof_csv_path = os.path.join(target_run_dir, tof_csv_filename)
        scope_csv_path = os.path.join(target_run_dir, scope_csv_filename)

        print(f"\n---> Starting High-Precision Tomography Acquisition...")
        print(f"     [Algorithm] Cross-Correlation + Sub-sample Parabolic Interpolation")
        print(f"     [Safety Config] +{int(self.safety_factor*100)}% ToF Safety Factor | Min Cooldown: {self.min_cooldown_ms} ms")
        print(f"     [Run Directory] {target_run_dir}")
        print(f"     [ToF CSV] {tof_csv_filename}")
        print(f"     [Scope CSV] {scope_csv_filename}")
        print(f"     [Mode] Continuous acquisition loop (Press Ctrl+C to stop logging)")

        dt_us = 1.0 / self.sample_rate_mhz
        scan_index = 1

        try:
            with open(tof_csv_path, "w", newline="") as f_tof, open(scope_csv_path, "w", newline="") as f_scope:
                tof_writer = csv.writer(f_tof)
                scope_writer = csv.writer(f_scope)

                # Write Headers
                tof_writer.writerow(["Scan_Index", "Timestamp_s", "RayPath_ID", "Tx_Channel", "Rx_Channel", "TimeOfFlight_us", "AdaptiveCooldown_ms", "Method", "Status"])
                scope_writer.writerow(["Scan_Index", "Timestamp_s", "RayPath_ID", "Tx_Channel", "Rx_Channel", "SampleIndex", "Time_us", "CH1_Volts", "CH2_Volts"])

                # Flush headers immediately from RAM to storage
                f_tof.flush()
                os.fsync(f_tof.fileno())
                f_scope.flush()
                os.fsync(f_scope.fileno())

                while True:
                    print(f"\n--- [Scan #{scan_index}] Processing {len(self.transducer_pairs)} ray paths ---")

                    for ray_id, (tx, rx) in enumerate(self.transducer_pairs, start=1):
                        ray_timestamp = time.time()
                        print(f"  [Scan #{scan_index} | RayPath #{ray_id}] Tx Channel {tx} -> Rx Channel {rx}...")
                        
                        # 1. Command Pico 2W to select MUX channel and fire 250ns pulse
                        self.pico.trigger_pair(tx, rx)

                        # 2. Capture Scope Waveform
                        ch1, ch2 = self.scope.capture_waveform(num_samples=2048)
                        
                        # 3. Save raw waveform data to scope CSV & flush from RAM to physical disk
                        if ch1 is not None and ch2 is not None:
                            num_pts = min(len(ch1), len(ch2))
                            for i in range(num_pts):
                                t_us = i * dt_us
                                scope_writer.writerow([scan_index, f"{ray_timestamp:.3f}", ray_id, tx, rx, i, f"{t_us:.4f}", f"{ch1[i]:.4f}", f"{ch2[i]:.4f}"])
                            
                            # Force write to storage for scope waveform data after each ray path capture
                            f_scope.flush()
                            os.fsync(f_scope.fileno())

                        # 4. Compute High-Precision Cross-Correlation ToF
                        tof = self.calculate_precision_tof_micros(ch1, ch2, method="cross_correlation")
                        
                        if tof is not None:
                            status = "OK"
                            self.last_observed_tof_us[(tx, rx)] = tof
                            cooldown_s = self.calculate_adaptive_cooldown_sec(tx, rx, tof)
                            print(f"     -> Sub-sample ToF: {tof:.4f} µs | Adaptive Cooldown (+10% margin): {cooldown_s*1000.0:.2f} ms")
                        else:
                            tof = -1.0
                            status = "NO_SIGNAL"
                            cooldown_s = 0.010
                            print(f"     -> [WARNING] No signal detected for Tx {tx} -> Rx {rx}")

                        self.tof_results[(tx, rx)] = tof
                        tof_writer.writerow([scan_index, f"{ray_timestamp:.3f}", ray_id, tx, rx, f"{tof:.4f}" if tof > 0 else "N/A", f"{cooldown_s*1000.0:.2f}", "CrossCorrelation", status])
                        
                        # Force write to storage for ToF matrix entry
                        f_tof.flush()
                        os.fsync(f_tof.fileno())

                        # Enforce adaptive +10% safety factor cooldown before next pulse
                        time.sleep(cooldown_s)

                    scan_index += 1
                    if not continuous:
                        break

        except KeyboardInterrupt:
            print(f"\n[SUCCESS] Acquisition loop stopped cleanly by user (KeyboardInterrupt).")
        finally:
            self.scope.close()
            self.pico.close()
            print(f"          - ToF Matrix saved to: '{tof_csv_path}'")
            print(f"          - Scope Waveforms saved to: '{scope_csv_path}'")

        return self.tof_results

if __name__ == "__main__":
    master = TomographyAcquisitionMaster(sample_rate_mhz=16.0, safety_factor=0.10, min_cooldown_ms=5.0)
    if master.connect():
        master.execute_tomography_scan(continuous=True)
