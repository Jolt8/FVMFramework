"""
AD9226 ADC Streaming Time-of-Flight Measurement
================================================

Replaces hantek_streaming_tof.py from the old Hantek 6022BE pipeline.

Since the Pico 2W now handles trigger synchronization in hardware
(PIO fires the TX pulse and starts ADC capture deterministically),
this script is MUCH simpler than the Hantek version:

  1. Command Pico to trigger + capture (single serial command)
  2. Receive the deterministically-timed ADC buffer
  3. Apply Hilbert envelope detection to find the echo peak
  4. Compute sub-sample ToF via parabolic interpolation
  5. Repeat for multiple loops, compute aggregate statistics

No more random USB timing, no more hoping the pulse lands in the
capture window, no more continuous streaming + software triggering.

Usage:
    python adc_streaming_tof.py
    python adc_streaming_tof.py --loops 10 --save-plot
    python adc_streaming_tof.py --min-tof 30 --max-tof 120
"""

import sys
import os
import time
import argparse
import numpy as np
from scipy.signal import hilbert

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pico_adc_interface import PicoADCInterface


class ADCStreamingToF:
    """
    Time-of-Flight measurement using the Pico 2W + AD9226 ADC.
    
    Key difference from the Hantek version: the trigger-to-sample timing
    is deterministic in hardware. The Pico starts the ADC clock, then
    fires the TX pulse a few microseconds later. The exact pre-trigger
    gap is small and consistent (MicroPython overhead for sm0.put(1)).
    
    We detect the TX pulse's effect on the waveform (or use the known
    pre-trigger offset) and measure the echo arrival relative to it.
    """
    
    def __init__(self, vref=1.0, pre_trigger_offset_us=5.0):
        """
        Args:
            vref: AD9226 reference voltage (ADC spans ±vref).
            pre_trigger_offset_us: Estimated time from ADC capture start
                                   to TX pulse firing. The Pico starts SM1,
                                   arms DMA, then fires SM0 — this takes
                                   ~5µs of MicroPython overhead.
        """
        self.adc = PicoADCInterface(vref=vref)
        self.pre_trigger_offset_us = pre_trigger_offset_us
    
    def connect(self):
        """Connects to the Pico 2W."""
        return self.adc.connect()
    
    def extract_tof(self, volts, dt_us, min_tof_us=20.0, max_tof_us=120.0,
                    noise_floor_factor=3.0):
        """
        Extracts Time-of-Flight from a single ADC capture using Hilbert envelope.
        
        Pipeline:
          1. Remove DC offset
          2. Estimate pre-trigger baseline noise level
          3. Compute Hilbert analytic signal envelope
          4. Search for envelope peak within [min_tof_us, max_tof_us] window
          5. Refine with parabolic sub-sample interpolation
          6. Subtract pre-trigger offset to get true ToF
        
        Args:
            volts: numpy array of voltage samples
            dt_us: time step in microseconds
            min_tof_us: minimum expected ToF (start of search window)
            max_tof_us: maximum expected ToF (end of search window)
            noise_floor_factor: envelope peak must exceed noise_floor * this factor
        
        Returns:
            (tof_us, peak_voltage, peak_time_us, snr) or (None, None, None, None)
        """
        n = len(volts)
        
        # 1. Remove DC offset
        v_ac = volts - np.mean(volts)
        
        # 2. Compute Hilbert envelope
        analytic = hilbert(v_ac)
        envelope = np.abs(analytic)
        
        # 3. Estimate noise floor from the first few samples (pre-trigger baseline)
        #    The first ~50 samples (5µs at 10 MSPS) are before the TX pulse
        baseline_end = min(50, n // 4)
        noise_floor = np.std(envelope[:baseline_end])
        
        # 4. Define search window
        #    The pre-trigger offset means the TX pulse fires at ~pre_trigger_offset_us
        #    into the capture. So the echo arrives at pre_trigger_offset_us + ToF.
        #    Search for the echo between:
        #      start = pre_trigger_offset_us + min_tof_us
        #      end   = pre_trigger_offset_us + max_tof_us
        search_start_us = self.pre_trigger_offset_us + min_tof_us
        search_end_us = self.pre_trigger_offset_us + max_tof_us
        
        search_start_idx = int(search_start_us / dt_us)
        search_end_idx = min(int(search_end_us / dt_us), n)
        
        if search_start_idx >= search_end_idx or search_start_idx >= n:
            return None, None, None, None
        
        # 5. Find envelope peak in search window
        search_region = envelope[search_start_idx:search_end_idx]
        if len(search_region) == 0:
            return None, None, None, None
        
        local_peak_idx = np.argmax(search_region)
        global_peak_idx = search_start_idx + local_peak_idx
        peak_envelope_v = envelope[global_peak_idx]
        
        # 6. SNR check — peak must be significantly above noise
        snr = peak_envelope_v / noise_floor if noise_floor > 1e-12 else float('inf')
        if snr < noise_floor_factor:
            return None, None, None, None
        
        # 7. Sub-sample parabolic interpolation on envelope
        delta = 0.0
        if 0 < local_peak_idx < len(search_region) - 1:
            y0 = search_region[local_peak_idx - 1]
            y1 = search_region[local_peak_idx]
            y2 = search_region[local_peak_idx + 1]
            denom = y0 - 2.0 * y1 + y2
            if abs(denom) > 1e-15:
                delta = 0.5 * (y0 - y2) / denom
        
        # 8. Calculate ToF
        #    Peak time in capture = (global_peak_idx + delta) * dt_us
        #    ToF = peak_time - pre_trigger_offset
        peak_time_us = (global_peak_idx + delta) * dt_us
        tof_us = peak_time_us - self.pre_trigger_offset_us
        
        peak_voltage = volts[global_peak_idx]
        
        if tof_us < 0:
            return None, None, None, None
        
        return tof_us, peak_voltage, peak_time_us, snr
    
    def measure(self, tx=1, rx=2, min_tof_us=20.0, max_tof_us=120.0, verbose=True):
        """
        Performs a single trigger-capture-extract cycle.
        
        Returns: dict with capture data and ToF result, or None on failure.
        """
        volts, sample_rate, dt_us = self.adc.trigger_and_capture(tx=tx, rx=rx)
        
        if volts is None:
            if verbose:
                print("[X] Capture failed")
            return None
        
        tof, peak_v, peak_t, snr = self.extract_tof(
            volts, dt_us, min_tof_us=min_tof_us, max_tof_us=max_tof_us
        )
        
        if verbose:
            if tof is not None:
                print(f"  ToF = {tof:.3f} µs | Echo peak = {peak_v:+.4f} V "
                      f"@ {peak_t:.2f} µs | SNR = {snr:.1f}")
            else:
                print(f"  [!] No valid echo detected (SNR too low or outside window)")
        
        return {
            'volts': volts,
            'sample_rate': sample_rate,
            'dt_us': dt_us,
            'tof_us': tof,
            'peak_voltage': peak_v,
            'peak_time_us': peak_t,
            'snr': snr,
            'tx': tx,
            'rx': rx,
        }
    
    def close(self):
        """Closes the serial connection."""
        self.adc.close()


def save_debug_plot(result, save_path):
    """Saves a 3-panel debug plot for visual verification."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    volts = result['volts']
    dt_us = result['dt_us']
    n = len(volts)
    t_us = np.arange(n) * dt_us
    
    v_ac = volts - np.mean(volts)
    envelope = np.abs(hilbert(v_ac))
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    fig.suptitle(f'AD9226 Streaming ToF — Tx {result["tx"]} → Rx {result["rx"]} '
                 f'@ {result["sample_rate"]/1e6:.1f} MSPS',
                 fontsize=14, fontweight='bold', y=0.98)
    
    # Panel 1: Raw waveform
    ax1 = axes[0]
    ax1.plot(t_us, volts, color='#1E90FF', linewidth=0.6, label='Raw ADC')
    if result['peak_time_us'] is not None:
        ax1.axvline(result['peak_time_us'], color='lime', linewidth=1.5,
                    linestyle='-.', label=f'Echo peak @ {result["peak_time_us"]:.2f} µs')
    ax1.set_title(f'Raw Waveform ({n} samples)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Voltage (V)')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Hilbert envelope
    ax2 = axes[1]
    ax2.plot(t_us, v_ac, color='#1E90FF', linewidth=0.5, alpha=0.4, label='AC Signal')
    ax2.plot(t_us, envelope, color='#FF4500', linewidth=1.5, label='Hilbert Envelope')
    ax2.plot(t_us, -envelope, color='#FF4500', linewidth=1.5, alpha=0.3)
    ax2.set_title('Hilbert Envelope Detection', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Voltage (V)')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: ToF measurement visualization
    ax3 = axes[2]
    ax3.plot(t_us, envelope, color='#FF4500', linewidth=1.2, label='Envelope')
    
    if result['tof_us'] is not None:
        ax3.axvline(result['peak_time_us'], color='lime', linewidth=2,
                    linestyle='-.', label=f'Echo peak: {result["peak_time_us"]:.2f} µs')
        ax3.set_title(f'ToF Measurement: {result["tof_us"]:.3f} µs | '
                      f'SNR: {result["snr"]:.1f}',
                      fontsize=11, fontweight='bold')
    else:
        ax3.set_title('No valid echo detected', fontsize=11, fontweight='bold')
    
    ax3.set_xlabel('Time (µs)')
    ax3.set_ylabel('Amplitude (V)')
    ax3.legend(loc='upper right')
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[+] Debug plot saved: {save_path}")


def main():
    parser = argparse.ArgumentParser(
        description='AD9226 ADC Streaming ToF Measurement',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python adc_streaming_tof.py
  python adc_streaming_tof.py --loops 10 --save-plot
  python adc_streaming_tof.py --min-tof 30 --max-tof 100
  python adc_streaming_tof.py --pre-trigger 3.0 --vref 0.5
        """
    )
    parser.add_argument('--loops', type=int, default=5,
                        help='Number of measurement loops (default: 5)')
    parser.add_argument('--min-tof', type=float, default=20.0,
                        help='Minimum ToF search window in µs (default: 20.0)')
    parser.add_argument('--max-tof', type=float, default=120.0,
                        help='Maximum ToF search window in µs (default: 120.0)')
    parser.add_argument('--pre-trigger', type=float, default=5.0,
                        help='Pre-trigger offset in µs (default: 5.0)')
    parser.add_argument('--vref', type=float, default=1.0,
                        help='AD9226 reference voltage in V (default: 1.0)')
    parser.add_argument('--tx', type=int, default=1,
                        help='TX transducer channel (default: 1)')
    parser.add_argument('--rx', type=int, default=2,
                        help='RX transducer channel (default: 2)')
    parser.add_argument('--save-plot', action='store_true',
                        help='Save a debug plot for the first measurement')
    parser.add_argument('--delay', type=float, default=0.05,
                        help='Delay between measurements in seconds (default: 0.05)')
    args = parser.parse_args()
    
    print("=" * 60)
    print("  AD9226 ADC STREAMING TIME-OF-FLIGHT MEASUREMENT")
    print("  (Hardware-deterministic PIO + DMA capture)")
    print("=" * 60)
    
    streamer = ADCStreamingToF(vref=args.vref, pre_trigger_offset_us=args.pre_trigger)
    
    if not streamer.connect():
        sys.exit(1)
    
    all_tof = []
    first_result = None
    
    try:
        for loop_idx in range(args.loops):
            print(f"\n{'=' * 60}")
            print(f"  Measurement {loop_idx + 1} / {args.loops}")
            print(f"{'=' * 60}")
            
            result = streamer.measure(
                tx=args.tx, rx=args.rx,
                min_tof_us=args.min_tof, max_tof_us=args.max_tof
            )
            
            if result and result['tof_us'] is not None:
                all_tof.append(result['tof_us'])
                if first_result is None:
                    first_result = result
            
            time.sleep(args.delay)
    
    except KeyboardInterrupt:
        print("\n[Stopped by user]")
    
    finally:
        streamer.close()
    
    # Aggregate statistics
    if len(all_tof) > 0:
        all_tof = np.array(all_tof)
        print(f"\n{'#' * 60}")
        print(f"  AGGREGATE RESULTS ({len(all_tof)} valid measurements)")
        print(f"  Mean ToF  : {np.mean(all_tof):.3f} µs")
        print(f"  Median ToF: {np.median(all_tof):.3f} µs")
        print(f"  Std Dev   : {np.std(all_tof):.3f} µs")
        print(f"  Min / Max : {np.min(all_tof):.3f} / {np.max(all_tof):.3f} µs")
        print(f"{'#' * 60}")
    else:
        print(f"\n[!] No valid ToF measurements obtained.")
        print("    Troubleshooting:")
        print("    1. Is the Pico running pico_adc_tomography.py?")
        print("    2. Are the transducers connected and in water?")
        print(f"    3. Try widening --min-tof / --max-tof (currently {args.min_tof}-{args.max_tof} µs)")
        print(f"    4. Try adjusting --pre-trigger (currently {args.pre_trigger} µs)")
        print(f"    5. Try adjusting --vref (currently {args.vref} V)")
    
    # Save debug plot
    if args.save_plot and first_result is not None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plot_path = os.path.join(script_dir, "streaming_tof_debug.png")
        try:
            save_debug_plot(first_result, plot_path)
        except ImportError:
            print("[*] matplotlib not available, skipping debug plot")


if __name__ == '__main__':
    main()
