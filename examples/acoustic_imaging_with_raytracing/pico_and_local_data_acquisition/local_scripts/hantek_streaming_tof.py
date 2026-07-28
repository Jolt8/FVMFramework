"""
hantek_streaming_tof.py — Continuous Streaming Software-Trigger ToF for Hantek 6022BE

This script replicates what OpenHantek does internally to achieve accurate,
consistent Time-of-Flight measurements:

  1. Streams continuous USB bulk data from the Hantek 6022BE
  2. Accumulates thousands of samples into a single contiguous numpy buffer
  3. Applies software edge-detection on CH1 to locate trigger pulse rising edges
  4. For each trigger, searches CH2 for the acoustic echo via Hilbert envelope detection
  5. Computes sub-sample ToF via parabolic interpolation on the envelope peak
  6. Averages across multiple pulse cycles for high-precision results

Why this works when capture_waveform() doesn't:
  The Hantek 6022BE has NO hardware trigger — the Cypress FX2 chip just continuously
  streams ADC samples into a USB FIFO. OpenHantek achieves its stable display by
  reading this stream continuously and doing software triggering on the CPU side.
  The old approach read a single 2048-byte frame and hoped the pulse landed in it,
  which is essentially random timing. This approach captures everything and finds
  the pulses in post-processing.

Prerequisites:
  - Hantek 6022BE connected via USB 2.0 port with WinUSB driver (via Zadig)
  - OpenHantek GUI must be CLOSED (exclusive USB access)
  - Pico 2W running better_pico_tomography.py in standalone continuous-fire mode

Usage:
  python hantek_streaming_tof.py
  python hantek_streaming_tof.py --duration 500 --loops 5 --save-plot
  python hantek_streaming_tof.py --rate-code 0x0003 --sample-rate 8.0
"""

import sys
import time
import os
import argparse
import numpy as np
from scipy.signal import hilbert
import libusb_package
import usb.core
import usb.util


class HantekStreamingToF:
    """
    Continuous streaming capture + software trigger for Hantek 6022BE.

    Unlike the single-frame capture_waveform() approach, this class:
    - Reads the USB bulk endpoint in a tight loop for a specified duration
    - Accumulates ALL data into a single contiguous buffer
    - Applies software triggering (like OpenHantek does) to find pulse edges
    - Extracts ToF from each pulse cycle and averages for precision
    """

    VID = 0x04B5
    PID = 0x6022
    BULK_EP = 0x86  # Bulk IN endpoint

    def __init__(self, rate_code = 0x0007, sample_rate_mhz = 1.0,
                 ch1_range_v = 5.0, ch2_range_v = 0.5,
                 pre_trigger_offset_us = 2.0):
        """
        Args:
            rate_code: USB control transfer wValue for sample rate (0xE2 command).
                       0x0007 = ~1 MS/s on the stock 6022BE firmware.
                       Try 0x0003 for higher rates if desired.
            sample_rate_mhz: Expected sample rate in MHz (for time axis calculations).
            ch1_range_v: CH1 full-scale voltage range (±V). 5.0 = ±5V.
            ch2_range_v: CH2 full-scale voltage range (±V). 0.5 = ±500mV.
            pre_trigger_offset_us: The PIO scope trigger (GP11) rises this many µs
                                   before the actual TX pulse fires on GP10.
                                   Subtracted from raw measured ToF.
        """
        self.dev = None
        self.rate_code = rate_code
        self.sample_rate_mhz = sample_rate_mhz
        self.dt_us = 1.0 / sample_rate_mhz
        self.ch1_range_v = ch1_range_v
        self.ch2_range_v = ch2_range_v
        self.pre_trigger_offset_us = pre_trigger_offset_us

    def connect(self):
        """Connect to Hantek 6022BE and configure gain + sample rate."""
        backend = libusb_package.get_libusb1_backend()
        self.dev = usb.core.find(idVendor = self.VID, idProduct = self.PID, backend = backend)

        if self.dev is None:
            print("[X] Hantek 6022BE not found on USB bus.")
            print("    Check: USB 2.0 port, WinUSB driver (Zadig), OpenHantek closed.")
            return False

        try:
            self.dev.set_configuration()
            usb.util.claim_interface(self.dev, 0)
        except usb.core.USBError as e:
            if "Access denied" in str(e) or e.errno == 13:
                print("[X] ACCESS DENIED — close OpenHantek GUI first.")
            else:
                print(f"[X] USB error: {e}")
            return False

        try:
            self.dev.set_interface_altsetting(0, 0)
        except Exception:
            pass

        # Set gain: CH1 = 5V (code 0x0001), CH2 = 500mV (code 0x0001) or 5V (0x0000)
        ch2_gain_code = 0x0001 if self.ch2_range_v <= 0.5 else 0x0000
        self.dev.ctrl_transfer(0x40, 0xE3, 0x0001, ch2_gain_code, b"\x01")

        # Set sample rate
        self.dev.ctrl_transfer(0x40, 0xE2, self.rate_code, 0x0000, b"\x01")

        # Set continuous trigger mode (same as test_bulk_loop.py / OpenHantek startup)
        try:
            self.dev.ctrl_transfer(0x40, 0xE4, 0x0001, 0x0000, b"\x01")
        except Exception:
            pass  # Some firmware versions don't support 0xE4

        print(f"[+] Hantek 6022BE connected")
        print(f"    Rate code: 0x{self.rate_code:04X} ({self.sample_rate_mhz} MS/s assumed)")
        print(f"    CH1: +/-{self.ch1_range_v}V | CH2: +/-{self.ch2_range_v}V")
        return True

    def stream_bulk(self, duration_ms = 200, chunk_bytes = 16384):
        """
        Streams continuous bulk USB data for the specified duration.

        This is the key difference from capture_waveform(): instead of reading
        a single 2048-byte frame, we read continuously for hundreds of ms,
        accumulating a large buffer that contains many complete pulse cycles.

        Returns: (raw_bytearray, n_reads, n_errors, elapsed_seconds)
        """
        # Flush stale data from endpoint FIFO
        try:
            self.dev.clear_halt(self.BULK_EP)
        except Exception:
            pass

        # Start ADC capture
        self.dev.ctrl_transfer(0x40, 0xE0, 0x0001, 0x0000, b"\x01")

        raw = bytearray()
        t_start = time.perf_counter()
        deadline = t_start + (duration_ms / 1000.0)
        read_count = 0
        error_count = 0

        while time.perf_counter() < deadline:
            try:
                data = self.dev.read(self.BULK_EP, chunk_bytes, timeout = 100)
                if data and len(data) > 0:
                    raw.extend(data)
                    read_count += 1
            except usb.core.USBTimeoutError:
                error_count += 1
                continue
            except usb.core.USBError:
                error_count += 1
                if error_count > 50:
                    break

        # Stop ADC capture
        try:
            self.dev.ctrl_transfer(0x40, 0xE0, 0x0000, 0x0000, b"\x01")
        except Exception:
            pass

        elapsed = time.perf_counter() - t_start
        return raw, read_count, error_count, elapsed

    def deinterleave(self, raw_bytes):
        """Split interleaved raw bytes into CH1 and CH2 voltage arrays."""
        n_pairs = len(raw_bytes) // 2
        if n_pairs < 10:
            return None, None

        data = np.frombuffer(raw_bytes[:n_pairs * 2], dtype = np.uint8)
        ch1_raw = data[0::2]
        ch2_raw = data[1::2]

        ch1_v = (ch1_raw.astype(np.float64) - 128.0) / 128.0 * self.ch1_range_v
        ch2_v = (ch2_raw.astype(np.float64) - 128.0) / 128.0 * self.ch2_range_v

        return ch1_v, ch2_v

    def find_trigger_edges(self, ch1_v, threshold_v = 1.5, min_spacing_samples = 50):
        """
        Software trigger: detects CH1 rising edges crossing threshold.

        This is what OpenHantek does internally — it continuously streams data
        and uses CPU-side edge detection to stabilize the display.

        Args:
            ch1_v: CH1 voltage array
            threshold_v: trigger threshold voltage
            min_spacing_samples: minimum samples between consecutive triggers
                                 (prevents re-triggering on the same pulse)

        Returns: array of sample indices where trigger rising edges occur
        """
        above = ch1_v > threshold_v
        transitions = np.diff(above.astype(np.int8))
        rising_edges = np.where(transitions == 1)[0] + 1

        if len(rising_edges) == 0:
            return np.array([], dtype = int)

        # Enforce minimum spacing between consecutive triggers
        filtered = [rising_edges[0]]
        for edge in rising_edges[1:]:
            if edge - filtered[-1] >= min_spacing_samples:
                filtered.append(edge)

        return np.array(filtered, dtype = int)

    def extract_tof_single(self, ch2_v, trig_idx,
                           min_tof_us = 20.0, max_tof_us = 120.0):
        """
        Extract ToF for a single trigger event using Hilbert envelope detection.

        Pipeline:
        1. Defines search window: [trig_idx + min_offset, trig_idx + max_offset]
        2. Removes DC offset from the window
        3. Computes Hilbert envelope (analytic signal magnitude)
        4. Finds envelope peak — more robust than raw peak, avoids cycle skipping
        5. Refines with parabolic sub-sample interpolation

        Returns: (tof_us, peak_voltage, peak_sample_idx) or (None, None, None)
        """
        min_offset = int(min_tof_us * self.sample_rate_mhz)
        max_offset = int(max_tof_us * self.sample_rate_mhz)

        win_start = trig_idx + min_offset
        win_end = trig_idx + max_offset

        if win_end >= len(ch2_v) or win_start >= len(ch2_v):
            return None, None, None

        window = ch2_v[win_start:win_end].copy()

        if len(window) < 10:
            return None, None, None

        # Remove DC offset from the search window
        window -= np.mean(window)

        # Hilbert envelope detection — gives the smooth amplitude envelope
        # of the wavepacket, eliminating cycle-skipping issues that plague
        # raw peak detection on oscillatory signals
        analytic_signal = hilbert(window)
        envelope = np.abs(analytic_signal)

        # Find envelope peak
        peak_local = np.argmax(envelope)
        peak_global = win_start + peak_local

        # Sub-sample parabolic interpolation on the envelope
        delta = 0.0
        if 0 < peak_local < len(envelope) - 1:
            y0 = envelope[peak_local - 1]
            y1 = envelope[peak_local]
            y2 = envelope[peak_local + 1]
            denom = y0 - 2.0 * y1 + y2
            if abs(denom) > 1e-15:
                delta = 0.5 * (y0 - y2) / denom

        # ToF = (echo_peak - trigger_edge) * dt, minus PIO pre-trigger offset
        tof_samples = (peak_global + delta) - trig_idx
        tof_us = tof_samples * self.dt_us - self.pre_trigger_offset_us

        peak_voltage = ch2_v[peak_global]

        return tof_us, peak_voltage, peak_global

    def measure(self, duration_ms = 200, trigger_threshold_v = 1.5,
                min_tof_us = 20.0, max_tof_us = 120.0, verbose = True):
        """
        Full measurement pipeline:
        1. Stream bulk data for duration_ms
        2. Software trigger on CH1 rising edges
        3. Extract ToF for each trigger via Hilbert envelope
        4. Return aggregated statistics and raw data
        """
        # --- Stream ---
        raw, n_reads, n_errors, elapsed_s = self.stream_bulk(duration_ms = duration_ms)
        ch1, ch2 = self.deinterleave(raw)

        if ch1 is None:
            print(f"[X] Insufficient data: {len(raw)} bytes from {n_reads} reads")
            return None

        total_us = len(ch1) * self.dt_us

        if verbose:
            data_rate_kbps = len(raw) / elapsed_s / 1024.0 if elapsed_s > 0 else 0
            print(f"\n[Stream] {len(raw):,} bytes -> {len(ch1):,} sample pairs "
                  f"({total_us:.0f} us = {total_us / 1000:.1f} ms)")
            print(f"         {n_reads} USB reads, {n_errors} timeouts, "
                  f"{elapsed_s * 1000:.0f} ms wall time ({data_rate_kbps:.0f} KB/s)")

        # --- Software Trigger ---
        # Minimum spacing: at least 100 us worth of samples to avoid
        # re-triggering on the same pulse (Pico fires every ~250-350 us)
        min_spacing = max(int(100.0 * self.sample_rate_mhz), 20)
        edges = self.find_trigger_edges(ch1, threshold_v = trigger_threshold_v,
                                        min_spacing_samples = min_spacing)

        if verbose:
            print(f"[Trigger] {len(edges)} rising edges found in CH1 "
                  f"(threshold = {trigger_threshold_v}V)")

        if len(edges) == 0:
            print("[!] No trigger edges found!")
            print("    -> Is the Pico running and firing pulses?")
            print("    -> Is CH1 connected to the scope trigger output (GP11)?")
            print(f"    CH1 stats: min = {np.min(ch1):.3f}V, max = {np.max(ch1):.3f}V, "
                  f"mean = {np.mean(ch1):.3f}V, std = {np.std(ch1):.3f}V")
            return {
                'ch1': ch1, 'ch2': ch2, 'trigger_edges': edges,
                'tof_values': np.array([]), 'peak_indices': np.array([], dtype = int),
                'peak_voltages': np.array([]),
                'n_triggers': 0, 'n_valid': 0,
                'mean_us': None, 'median_us': None, 'std_us': None,
                'dt_us': self.dt_us,
            }

        # Diagnostic: trigger interval consistency (sanity check on sample rate)
        if len(edges) >= 2 and verbose:
            intervals_us = np.diff(edges) * self.dt_us
            print(f"         Trigger interval: {np.mean(intervals_us):.1f} +/- "
                  f"{np.std(intervals_us):.1f} us "
                  f"(expect ~250-350 us from Pico)")
            # If the measured interval is wildly off, the sample rate might be wrong
            if np.mean(intervals_us) < 50 or np.mean(intervals_us) > 5000:
                print(f"    [WARNING] Trigger interval seems off — double check "
                      f"--sample-rate ({self.sample_rate_mhz} MS/s) and "
                      f"--rate-code (0x{self.rate_code:04X})")

        # --- Extract ToF for each trigger ---
        tof_list = []
        peak_voltages = []
        peak_indices = []

        for i, edge in enumerate(edges):
            tof, peak_v, peak_idx = self.extract_tof_single(
                ch2, edge, min_tof_us = min_tof_us, max_tof_us = max_tof_us
            )
            if tof is not None and tof > 0:
                tof_list.append(tof)
                peak_voltages.append(peak_v)
                peak_indices.append(peak_idx)
                if verbose:
                    edge_t = edge * self.dt_us
                    print(f"  Edge #{i + 1:3d} @ {edge_t:8.1f} us -> "
                          f"ToF = {tof:7.3f} us  (echo {peak_v:+.4f} V)")

        tof_arr = np.array(tof_list)

        if verbose:
            print(f"\n{'=' * 60}")
            print(f"  Valid measurements: {len(tof_arr)} / {len(edges)} triggers")
            if len(tof_arr) > 0:
                print(f"  Mean ToF  : {np.mean(tof_arr):.3f} us")
                print(f"  Median ToF: {np.median(tof_arr):.3f} us")
                print(f"  Std Dev   : {np.std(tof_arr):.3f} us")
                print(f"  Min / Max : {np.min(tof_arr):.3f} / {np.max(tof_arr):.3f} us")
            else:
                print("  [!] No valid ToF extracted from any trigger")
            print(f"{'=' * 60}")

        return {
            'ch1': ch1,
            'ch2': ch2,
            'trigger_edges': edges,
            'tof_values': tof_arr,
            'peak_voltages': np.array(peak_voltages),
            'peak_indices': np.array(peak_indices, dtype = int) if peak_indices else np.array([], dtype = int),
            'n_triggers': len(edges),
            'n_valid': len(tof_arr),
            'mean_us': float(np.mean(tof_arr)) if len(tof_arr) > 0 else None,
            'median_us': float(np.median(tof_arr)) if len(tof_arr) > 0 else None,
            'std_us': float(np.std(tof_arr)) if len(tof_arr) > 0 else None,
            'dt_us': self.dt_us,
        }

    def close(self):
        """Release USB interface and stop capture."""
        if self.dev:
            try:
                self.dev.ctrl_transfer(0x40, 0xE0, 0x0000, 0x0000, b"\x01")
            except Exception:
                pass
            try:
                usb.util.release_interface(self.dev, 0)
            except Exception:
                pass
            self.dev = None
            print("[+] USB interface released.")


def save_debug_plot(result, dt_us, save_path, max_triggers_shown = 8):
    """Saves a 3-panel matplotlib debug plot for visual verification."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ch1 = result['ch1']
    ch2 = result['ch2']
    edges = result['trigger_edges']
    tof_arr = result['tof_values']
    peak_indices = result['peak_indices']

    t_us = np.arange(len(ch1)) * dt_us

    fig, axes = plt.subplots(3, 1, figsize = (14, 12))
    fig.suptitle('Hantek 6022BE Streaming ToF — Software Trigger Debug',
                 fontsize = 14, fontweight = 'bold', y = 0.98)

    # ---- Panel 1: Full stream overview ----
    ax1 = axes[0]
    ax1.plot(t_us, ch1, color = '#FFD700', linewidth = 0.4, alpha = 0.8, label = 'CH1 (Trigger)')
    ax1.plot(t_us, ch2, color = '#1E90FF', linewidth = 0.4, alpha = 0.8, label = 'CH2 (Echo)')
    for i, edge in enumerate(edges[:max_triggers_shown]):
        ax1.axvline(edge * dt_us, color = 'red', linewidth = 0.7, alpha = 0.5,
                    label = 'Trigger edges' if i == 0 else None)
    ax1.set_title(f'Full Stream: {len(ch1):,} samples, {len(edges)} triggers found',
                  fontsize = 11, fontweight = 'bold')
    ax1.set_ylabel('Voltage (V)')
    ax1.legend(loc = 'upper right', fontsize = 9)
    ax1.grid(True, alpha = 0.3)

    # ---- Panel 2: Zoom into first trigger+echo cycle ----
    ax2 = axes[1]
    if len(edges) > 0:
        first_edge = edges[0]
        zoom_start = max(0, first_edge - int(10.0 / dt_us))
        zoom_end = min(len(ch1), first_edge + int(140.0 / dt_us))
        zoom_t = np.arange(zoom_start, zoom_end) * dt_us

        ax2.plot(zoom_t, ch1[zoom_start:zoom_end], color = '#FFD700', linewidth = 1.2,
                 label = 'CH1 (Trigger)')
        ax2.plot(zoom_t, ch2[zoom_start:zoom_end], color = '#1E90FF', linewidth = 1.2,
                 label = 'CH2 (Echo)')
        ax2.axvline(first_edge * dt_us, color = 'red', linewidth = 1.5, linestyle = '--',
                    label = f'Trigger @ {first_edge * dt_us:.1f} us')

        if len(peak_indices) > 0 and len(tof_arr) > 0:
            ax2.axvline(peak_indices[0] * dt_us, color = 'lime', linewidth = 1.5,
                        linestyle = '-.',
                        label = f'Echo peak @ {peak_indices[0] * dt_us:.1f} us')
            ax2.set_title(f'First Trigger Zoom — Measured ToF = {tof_arr[0]:.3f} us',
                          fontsize = 11, fontweight = 'bold')
        else:
            ax2.set_title('First Trigger Zoom — No echo detected',
                          fontsize = 11, fontweight = 'bold')

        ax2.set_ylabel('Voltage (V)')
        ax2.legend(loc = 'upper right', fontsize = 9)
        ax2.grid(True, alpha = 0.3)
    else:
        ax2.text(0.5, 0.5, 'No triggers found', transform = ax2.transAxes,
                 ha = 'center', va = 'center', fontsize = 14)

    # ---- Panel 3: Hilbert envelope of the search window ----
    ax3 = axes[2]
    if len(edges) > 0 and len(tof_arr) > 0:
        first_edge = edges[0]
        min_offset = int(20.0 / dt_us)
        max_offset = int(120.0 / dt_us)
        win_start = first_edge + min_offset
        win_end = min(first_edge + max_offset, len(ch2))

        if win_end > win_start + 10 and win_end <= len(ch2):
            window = ch2[win_start:win_end].copy()
            window -= np.mean(window)
            envelope = np.abs(hilbert(window))
            win_t = np.arange(win_start, win_end) * dt_us

            ax3.plot(win_t, window, color = '#1E90FF', linewidth = 0.8, alpha = 0.5,
                     label = 'CH2 (DC removed)')
            ax3.plot(win_t, envelope, color = '#FF4500', linewidth = 1.8,
                     label = 'Hilbert Envelope')
            ax3.plot(win_t, -envelope, color = '#FF4500', linewidth = 1.8, alpha = 0.3)

            peak_local = np.argmax(envelope)
            peak_t = (win_start + peak_local) * dt_us
            ax3.axvline(peak_t, color = 'lime', linewidth = 1.5, linestyle = '-.',
                        label = f'Envelope peak @ {peak_t:.1f} us')

            ax3.set_title('Hilbert Envelope Detection (search window)',
                          fontsize = 11, fontweight = 'bold')

        ax3.set_ylabel('Voltage (V)')
        ax3.legend(loc = 'upper right', fontsize = 9)
        ax3.grid(True, alpha = 0.3)
    else:
        ax3.text(0.5, 0.5, 'No valid ToF data', transform = ax3.transAxes,
                 ha = 'center', va = 'center', fontsize = 14)

    axes[-1].set_xlabel('Time (us)')

    plt.tight_layout(rect = [0, 0, 1, 0.96])
    plt.savefig(save_path, dpi = 150, bbox_inches = 'tight')
    plt.close()
    print(f"[+] Debug plot saved: {save_path}")


# ---------------------------------------------------------------------------
# CLI Entry Point
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description = 'Hantek 6022BE Continuous Streaming ToF Measurement',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        epilog = """
Examples:
  python hantek_streaming_tof.py
  python hantek_streaming_tof.py --duration 500 --loops 5
  python hantek_streaming_tof.py --rate-code 0x0003 --sample-rate 8.0
  python hantek_streaming_tof.py --save-plot --min-tof 40 --max-tof 100
        """
    )
    parser.add_argument('--rate-code', type = lambda x: int(x, 0), default = 0x0007,
                        help = 'Sample rate code for 0xE2 command (default: 0x0007 = ~1 MS/s)')
    parser.add_argument('--sample-rate', type = float, default = 1.0,
                        help = 'Expected sample rate in MS/s (default: 1.0)')
    parser.add_argument('--duration', type = int, default = 200,
                        help = 'Stream duration per measurement in ms (default: 200)')
    parser.add_argument('--loops', type = int, default = 5,
                        help = 'Number of measurement loops (default: 5)')
    parser.add_argument('--threshold', type = float, default = 1.5,
                        help = 'CH1 trigger threshold in V (default: 1.5)')
    parser.add_argument('--min-tof', type = float, default = 20.0,
                        help = 'Minimum ToF search window in us (default: 20.0)')
    parser.add_argument('--max-tof', type = float, default = 120.0,
                        help = 'Maximum ToF search window in us (default: 120.0)')
    parser.add_argument('--pre-trigger-offset', type = float, default = 2.0,
                        help = 'PIO pre-trigger offset in us (default: 2.0)')
    parser.add_argument('--ch2-range', type = float, default = 0.5,
                        help = 'CH2 voltage range in V (default: 0.5 = +/-500mV)')
    parser.add_argument('--save-plot', action = 'store_true',
                        help = 'Save a matplotlib debug plot for the first measurement')
    args = parser.parse_args()

    print("=" * 60)
    print("  HANTEK 6022BE CONTINUOUS STREAMING ToF MEASUREMENT")
    print("  (Software Trigger — OpenHantek-style data pipeline)")
    print("=" * 60)

    scope = HantekStreamingToF(
        rate_code = args.rate_code,
        sample_rate_mhz = args.sample_rate,
        ch1_range_v = 5.0,
        ch2_range_v = args.ch2_range,
        pre_trigger_offset_us = args.pre_trigger_offset,
    )

    if not scope.connect():
        sys.exit(1)

    # Discard first capture to flush stale FIFO data from previous sessions
    print("\n[*] Flushing stale FIFO data...")
    scope.stream_bulk(duration_ms = 50)
    time.sleep(0.02)

    all_tof = []
    first_result = None

    try:
        for loop_idx in range(args.loops):
            print(f"\n{'=' * 60}")
            print(f"  Measurement {loop_idx + 1} / {args.loops}")
            print(f"{'=' * 60}")

            result = scope.measure(
                duration_ms = args.duration,
                trigger_threshold_v = args.threshold,
                min_tof_us = args.min_tof,
                max_tof_us = args.max_tof,
            )

            if result and len(result['tof_values']) > 0:
                all_tof.extend(result['tof_values'])
                if first_result is None:
                    first_result = result

            # Small gap between loops to let USB settle
            time.sleep(0.05)

    except KeyboardInterrupt:
        print("\n[Stopped by user]")

    finally:
        scope.close()

    # --- Final Aggregate Statistics ---
    if len(all_tof) > 0:
        all_tof = np.array(all_tof)
        print(f"\n{'#' * 60}")
        print(f"  AGGREGATE RESULTS ({len(all_tof)} total measurements)")
        print(f"  Mean ToF  : {np.mean(all_tof):.3f} us")
        print(f"  Median ToF: {np.median(all_tof):.3f} us")
        print(f"  Std Dev   : {np.std(all_tof):.3f} us")
        print(f"  Min / Max : {np.min(all_tof):.3f} / {np.max(all_tof):.3f} us")
        print(f"{'#' * 60}")
    else:
        print("\n[!] No valid ToF measurements obtained across all loops.")
        print("    Troubleshooting:")
        print("    1. Is the Pico running better_pico_tomography.py (standalone fire mode)?")
        print("    2. Is CH1 connected to GP11 (scope trigger output)?")
        print("    3. Is CH2 connected to the receive transducer?")
        print(f"    4. Try lowering --threshold (currently {args.threshold}V)")
        print("    5. Try widening --min-tof / --max-tof search window")

    # --- Save Debug Plot ---
    if args.save_plot and first_result is not None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plot_path = os.path.join(script_dir, "streaming_tof_debug.png")
        save_debug_plot(first_result, scope.dt_us, plot_path)


if __name__ == '__main__':
    main()
