"""
AD9226 ADC Capture Sanity Test
==============================

Connects to the Pico 2W, fires a single TX pulse, captures the ADC response,
prints waveform statistics, and saves a debug plot.

Equivalent to test_sync_order.py from the old Hantek pipeline, but much simpler
since the Pico now handles trigger synchronization in hardware.

Usage:
    python test_adc_capture.py
"""

import os
import sys
import numpy as np

# Add the local_scripts directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pico_adc_interface import PicoADCInterface


def save_debug_plot(volts, sample_rate, dt_us, save_path, tx=1, rx=2):
    """Saves a 3-panel debug plot of the captured waveform."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.signal import hilbert

    n = len(volts)
    t_us = np.arange(n) * dt_us

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 11))
    fig.suptitle(f'AD9226 ADC Capture — Tx {tx} → Rx {rx} @ {sample_rate/1e6:.1f} MSPS',
                 fontsize=14, fontweight='bold', y=0.98)

    # --- Panel 1: Raw Waveform ---
    ax1.plot(t_us, volts, color='#1E90FF', linewidth=0.8, label='ADC Voltage')
    ax1.set_title(f'Raw ADC Waveform ({n} samples, {t_us[-1]:.1f} µs window)',
                  fontsize=11, fontweight='bold')
    ax1.set_ylabel('Voltage (V)')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')

    # --- Panel 2: DC-Removed + Hilbert Envelope ---
    v_ac = volts - np.mean(volts)
    try:
        envelope = np.abs(hilbert(v_ac))
    except Exception:
        envelope = np.abs(v_ac)

    ax2.plot(t_us, v_ac, color='#1E90FF', linewidth=0.6, alpha=0.5, label='AC Signal')
    ax2.plot(t_us, envelope, color='#FF4500', linewidth=1.5, label='Hilbert Envelope')
    ax2.plot(t_us, -envelope, color='#FF4500', linewidth=1.5, alpha=0.3)
    ax2.set_title('DC-Removed Signal + Hilbert Envelope', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Voltage (V)')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')

    # --- Panel 3: Amplitude (Envelope Peak Detection) ---
    ax3.plot(t_us, envelope, color='#FF4500', linewidth=1.2, label='Envelope Amplitude')
    
    # Find and mark the peak
    peak_idx = np.argmax(envelope)
    peak_t = t_us[peak_idx]
    peak_v = envelope[peak_idx]
    ax3.axvline(peak_t, color='lime', linewidth=1.5, linestyle='-.',
                label=f'Peak @ {peak_t:.2f} µs ({peak_v:.4f} V)')
    ax3.scatter([peak_t], [peak_v], color='lime', s=80, zorder=5)
    
    ax3.set_title(f'Envelope Peak Detection — Peak at {peak_t:.2f} µs',
                  fontsize=11, fontweight='bold')
    ax3.set_xlabel('Time (µs)')
    ax3.set_ylabel('Amplitude (V)')
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc='upper right')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[+] Debug plot saved: {save_path}")


def main():
    print("=" * 60)
    print("  AD9226 ADC CAPTURE SANITY TEST")
    print("=" * 60)

    adc = PicoADCInterface()
    if not adc.connect():
        print("[X] Could not connect to Pico 2W.")
        return

    # Get status
    status = adc.get_status()
    if status:
        print(f"[+] Pico status: {status}")

    # Capture
    tx, rx = 1, 2
    print(f"\n[*] Triggering Tx {tx} -> Rx {rx}...")
    volts, sample_rate, dt_us = adc.trigger_and_capture(tx=tx, rx=rx)

    if volts is None:
        print("[X] Capture failed!")
        adc.close()
        return

    # Statistics
    n = len(volts)
    t_total_us = n * dt_us
    
    print(f"\n{'=' * 60}")
    print(f"  CAPTURE RESULTS")
    print(f"{'=' * 60}")
    print(f"  Samples     : {n}")
    print(f"  Sample Rate  : {sample_rate/1e6:.1f} MSPS")
    print(f"  Time Step    : {dt_us:.3f} µs ({1/dt_us:.0f} ns)")
    print(f"  Window       : {t_total_us:.1f} µs")
    print(f"  V min        : {np.min(volts):.4f} V")
    print(f"  V max        : {np.max(volts):.4f} V")
    print(f"  V mean       : {np.mean(volts):.4f} V")
    print(f"  V std        : {np.std(volts):.4f} V")
    print(f"  V pk-pk      : {np.ptp(volts):.4f} V")
    
    # Check for clipping (ADC at rail)
    n_clipped_low = np.sum(volts <= -0.99 * adc.vref)
    n_clipped_high = np.sum(volts >= 0.99 * adc.vref)
    if n_clipped_low > 0 or n_clipped_high > 0:
        print(f"  [WARNING] Clipping detected: {n_clipped_low} low, {n_clipped_high} high")
    
    # Check for signal presence
    v_ac = volts - np.mean(volts)
    if np.std(v_ac) < 0.001:
        print(f"  [WARNING] Very low signal amplitude (std < 1 mV). Check connections.")
    else:
        # Find dominant frequency via FFT
        fft_mag = np.abs(np.fft.rfft(v_ac))
        freqs = np.fft.rfftfreq(n, d=dt_us * 1e-6)
        if len(fft_mag) > 1:
            # Ignore DC bin
            peak_freq_idx = np.argmax(fft_mag[1:]) + 1
            peak_freq = freqs[peak_freq_idx]
            print(f"  Dominant Freq: {peak_freq/1e6:.3f} MHz")
    
    print(f"{'=' * 60}")

    # Save debug plot
    script_dir = os.path.dirname(os.path.abspath(__file__))
    plot_path = os.path.join(script_dir, "adc_capture_debug.png")
    
    try:
        save_debug_plot(volts, sample_rate, dt_us, plot_path, tx=tx, rx=rx)
    except ImportError:
        print("[*] matplotlib/scipy not available, skipping debug plot")

    adc.close()
    print("\n[+] Test complete.")


if __name__ == "__main__":
    main()
