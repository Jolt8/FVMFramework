"""
Hantek 6022BE Python Data Acquisition & Ultrasound Tomography Interface
Direct USB communication via PyUSB & libusb. Saves data to CSV.

TROUBLESHOOTING & RECOVERY NOTES:
---------------------------------
1. Device not showing up in Zadig ("List All Devices" checked, but missing):
   - Cause: Hardware connection dropped / USB 3.0 controller handshake issue with Cypress FX2 chip
     (Windows PnP status shows Present: False / CM_PROB_PHANTOM).
   - Fix: Plug directly into a USB 2.0 port on the motherboard (avoid USB 3.0 hubs or extension cables).

2. OpenHantek crashing on launch (ucrtbase.dll / Exception 0xc0000409 / BEX64):
   - Cause: Corrupted window geometry or scope state saved in Windows Registry under HKCU\Software\OpenHantek.
   - Fix: Run PowerShell command to wipe the corrupted registry key:
     Remove-Item -Path "HKCU:\Software\OpenHantek" -Recurse -Force
   - Ensure target driver in Zadig is set to WinUSB (v6.1+) for USB\VID_04B5&PID_6022.
"""


import sys
import time
import csv
import numpy as np
import libusb_package
import usb.core
import usb.util

class Hantek6022BE:
    VID_FW_LOADED  = 0x04B5
    PID_FW_LOADED  = 0x6022
    VID_BOOTLOADER = 0x04B4
    PID_BOOTLOADER = 0x6022

    GAIN_MAP = {
        5.0: 0x01,  # +/- 5V range
        2.0: 0x02,  # +/- 2V range
        1.0: 0x05,  # +/- 1V range
        0.5: 0x0A,  # +/- 500mV range
    }

    def __init__(self):
        self.dev = None
        self.ch1_gain_v = 5.0
        self.ch2_gain_v = 5.0
        self.sample_rate_mhz = 16.0

    def find_device(self):
        """Scans USB bus for Hantek 6022BE device using LibUSB backend."""
        try:
            backend = libusb_package.get_libusb1_backend()
            return usb.core.find(idVendor=self.VID_FW_LOADED, idProduct=self.PID_FW_LOADED, backend=backend)
        except Exception:
            for dev in libusb_package.find(find_all=True):
                if dev.idVendor in (self.VID_FW_LOADED, self.VID_BOOTLOADER) and dev.idProduct == self.PID_FW_LOADED:
                    return dev
        return None

    def connect(self):
        """Connects to the Hantek 6022BE USB scope."""
        self.dev = self.find_device()
        if self.dev is None:
            print("[X] Hantek 6022BE device not found on USB bus.")
            return False

        print(f"[+] Found Hantek 6022BE (VID=0x{self.dev.idVendor:04X}, PID=0x{self.dev.idProduct:04X})")

        try:
            self.dev.set_configuration()
            usb.util.claim_interface(self.dev, 0)
            print("[+] Claimed USB Interface #0")
            return True
        except usb.core.USBError as e:
            if "Access denied" in str(e) or e.errno == 13:
                print("\n[!] ACCESS DENIED: Another process (like OpenHantek GUI) is holding the USB handle.")
                print("    Please close OpenHantek GUI to allow direct Python USB access.\n")
            else:
                print(f"[X] USB Connection Error: {e}")
            return False

    def configure(self, sample_rate_mhz=16.0, ch1_range_v=5.0, ch2_range_v=5.0):
        """Configures sample rate and input gain."""
        if not self.dev:
            return False

        self.ch1_gain_v = ch1_range_v
        self.ch2_gain_v = ch2_range_v
        self.sample_rate_mhz = sample_rate_mhz

        # Send vendor commands to configure FPGA ADC clock and gain
        try:
            # Set Gain (Request 0xE3): CH1=1 (5V), CH2=1 (5V)
            self.dev.ctrl_transfer(0x40, 0xE3, 0x0001, 0x0001, b"\x01")
            # Set Sample Rate (Request 0xE2)
            self.dev.ctrl_transfer(0x40, 0xE2, 0x0001, 0x0000, b"\x01")
        except Exception:
            pass

        return True

    def capture_waveform(self, num_samples=2048):
        """Captures a clean 244 µs contiguous single-frame scope capture at 1 MS/s."""
        if not self.dev:
            return None, None

        try:
            self.dev.set_interface_altsetting(0, 0)
            time.sleep(0.005)
            
            ch2_val = 0x0001 if self.ch2_gain_v <= 0.5 else 0x0000
            # Set Gain (0xE3): CH1=5V (1), CH2=500mV (1) or 5V (0)
            self.dev.ctrl_transfer(0x40, 0xE3, 0x0001, ch2_val, b"\x01")
            # Set Samplerate (0xE2): Code 7 = 1 MS/s (244µs window)
            self.dev.ctrl_transfer(0x40, 0xE2, 0x0007, 0x0000, b"\x01")
            # Start hardware capture (0xE0)
            self.dev.ctrl_transfer(0x40, 0xE0, 0x0001, 0x0000, b"\x01")
            
            self.dev.clear_halt(0x86)
            data = self.dev.read(0x86, 2048, timeout=200)
            if not data or len(data) < 100:
                return None, None

            total_samples = len(data) // 2
            data_np = np.frombuffer(data[:total_samples * 2], dtype=np.uint8)
            ch1_raw = data_np[0::2]
            ch2_raw = data_np[1::2]

            ch1_volts = (ch1_raw.astype(np.float32) - 128.0) / 128.0 * self.ch1_gain_v
            ch2_volts = (ch2_raw.astype(np.float32) - 128.0) / 128.0 * self.ch2_gain_v

            return ch1_volts, ch2_volts
        except Exception:
            return None, None

    def save_to_csv(self, ch1_volts, ch2_volts, filename="capture_hantek_main.csv"):
        """Saves captured channel waveforms to a CSV file."""
        if ch1_volts is None or len(ch1_volts) == 0:
            print("[X] Cannot save CSV: No waveform data.")
            return False

        dt_us = 1.0 / self.sample_rate_mhz
        num_pts = min(len(ch1_volts), len(ch2_volts))

        print(f"[+] Writing {num_pts} samples to '{filename}'...")
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['SampleIndex', 'Time_us', 'CH1_Volts', 'CH2_Volts'])
            for i in range(num_pts):
                t_us = i * dt_us
                writer.writerow([i, f"{t_us:.4f}", f"{ch1_volts[i]:.4f}", f"{ch2_volts[i]:.4f}"])

        print(f"[SUCCESS] Saved CSV file: '{filename}'")
        return True

    def close(self):
        if self.dev:
            try:
                usb.util.release_interface(self.dev, 0)
                print("[+] USB interface released.")
            except Exception:
                pass
            self.dev = None

def main():
    print("=== Hantek 6022BE Live Capture & CSV Export ===")
    scope = Hantek6022BE()
    if scope.connect():
        scope.configure(sample_rate_mhz=16.0, ch1_range_v=5.0, ch2_range_v=5.0)
        print("\nCapturing pulse waveform from scope...")
        ch1, ch2 = scope.capture_waveform(num_samples=2048)
        
        if ch1 is not None and len(ch1) > 0:
            scope.save_to_csv(ch1, ch2, "capture_hantek_main.csv")
            
            # Time of Arrival (ToA) calculation for tomography
            thresh = 0.2  # 200mV threshold
            crossings = np.where(np.abs(ch1) > thresh)[0]
            if len(crossings) > 0:
                toa_us = crossings[0] / scope.sample_rate_mhz
                print(f"[TOMOGRAPHY] Pulse arrival detected at sample #{crossings[0]} ({toa_us:.3f} µs)")
            else:
                print("[TOMOGRAPHY] Baseline signal (no pulse > 200mV detected).")

        scope.close()

if __name__ == "__main__":
    main()
