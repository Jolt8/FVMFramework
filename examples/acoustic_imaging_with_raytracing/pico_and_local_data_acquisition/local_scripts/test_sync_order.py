import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from tomography_master_acquisition import TomographyAcquisitionMaster

def test_correct_sync_order():
    artifact_dir = r"C:\Users\wille\.gemini\antigravity-ide\brain\f196533a-4d85-4b5f-b49a-8df8c9df96f9"
    os.makedirs(artifact_dir, exist_ok=True)
    plot_path = os.path.join(artifact_dir, "scope_algorithm_debug.png")

    print("=== Testing Correct Scope-First Synchronization Order ===")
    master = TomographyAcquisitionMaster()
    if not master.connect():
        print("[X] Failed to connect.")
        return

    # STEP 1: Issue 0xE0 conversion start command to Hantek scope FIRST
    try:
        master.scope.dev.set_interface_altsetting(0, 0)
        master.scope.dev.ctrl_transfer(0x40, 0xE3, 0x0001, 0x0001, b"\x01")
        master.scope.dev.ctrl_transfer(0x40, 0xE2, 0x0003, 0x0000, b"\x01")
        master.scope.dev.ctrl_transfer(0x40, 0xE0, 0x0001, 0x0000, b"\x01")
    except Exception as e:
        print(f"[*] Scope start warning: {e}")

    # STEP 2: Fire Pico pulse IMMEDIATELY while scope is recording!
    t0 = time.perf_counter()
    master.pico.trigger_pair(1, 2)
    t1 = time.perf_counter()

    # STEP 3: Read USB bulk buffer from scope
    raw_bytes = bytearray()
    num_bursts = 8
    for b in range(num_bursts):
        try:
            master.scope.dev.clear_halt(0x86)
            data = master.scope.dev.read(0x86, 2048, timeout=200)
            if data and len(data) > 0:
                raw_bytes.extend(data)
            master.scope.dev.ctrl_transfer(0x40, 0xE0, 0x0001, 0x0000, b"\x01")
        except Exception:
            pass

    master.scope.close()
    master.pico.close()

    total_samples = len(raw_bytes) // 2
    if total_samples < 200:
        print(f"[X] Insufficient raw bytes ({len(raw_bytes)} bytes)")
        return

    data_np = np.frombuffer(raw_bytes[:total_samples * 2], dtype=np.uint8)
    ch1_raw = data_np[0::2]
    ch2_raw = data_np[1::2]

    ch1 = (ch1_raw.astype(np.float32) - 128.0) / 128.0 * 5.0
    ch2 = (ch2_raw.astype(np.float32) - 128.0) / 128.0 * 0.5

    dt_us = 1.0 / 16.0
    t_us = np.arange(total_samples) * dt_us

    # 1. Search for CH1 trigger edge (> 1.5 V)
    ch1_trig = np.where(ch1 > 1.5)[0]
    if len(ch1_trig) > 0:
        trig_idx = ch1_trig[0]
        trig_t = t_us[trig_idx]
        trig_found = True
    else:
        trig_idx = 0
        trig_t = 0.0
        trig_found = False

    # 2. Moving baseline detrending on CH2
    kernel_size = 25
    baseline = np.convolve(ch2, np.ones(kernel_size)/kernel_size, mode='same')
    ch2_detrended = ch2 - baseline

    # 3. Peak Search
    min_search_samples = int(15.0 * 16.0)
    search_start = min(trig_idx + min_search_samples, total_samples - 1) if trig_found else 0
    search_region = ch2_detrended[search_start:]

    if len(search_region) > 0:
        local_peak_i = np.argmax(np.abs(search_region))
        global_peak_idx = search_start + local_peak_i
        peak_t = global_peak_idx * dt_us
        calculated_tof = (global_peak_idx - trig_idx) * dt_us
        peak_v = ch2[global_peak_idx]
    else:
        global_peak_idx = 0
        peak_t = 0.0
        calculated_tof = -1.0
        peak_v = 0.0

    print(f"\n=======================================================")
    print(f" Pico Trigger Roundtrip Time: {(t1-t0)*1000.0:.2f} ms")
    print(f" Captured {total_samples} samples ({t_us[-1]:.2f} µs total window)")
    print(f" CH1 Trigger Edge Found: {trig_found} at sample #{trig_idx} (t = {trig_t:.3f} µs)")
    print(f" CH2 Echo Peak Found at sample #{global_peak_idx} (t = {peak_t:.3f} µs)")
    print(f" Calculated ToF: {calculated_tof:.3f} µs | Echo Peak Volts: {peak_v:.3f} V")
    print(f"=======================================================\n")

    # Save 3-Panel Debug Plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(11, 10))

    ax1.plot(t_us, ch1, label="CH1 (TX Trigger Ref)", color="#1f77b4", linewidth=1.2)
    ax1.plot(t_us, ch2, label="CH2 (RX Acoustic Raw)", color="#ff7f0e", linewidth=1.2)
    if trig_found:
        ax1.axvline(trig_t, color="blue", linestyle="--", linewidth=1.5, label=f"CH1 Trigger Edge @ {trig_t:.2f} µs")
    ax1.set_title(f"Scope-First Sync Order Raw Waveforms ({total_samples} samples @ 16 MS/s)", fontsize=11, fontweight='bold')
    ax1.set_ylabel("Voltage (V)")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper right")

    ax2.plot(t_us, ch2_detrended, label="CH2 Detrended Signal (RC Ramp Filtered)", color="#2ca02c", linewidth=1.2)
    ax2.plot(t_us, baseline, label="Extracted RC Baseline", color="gray", linestyle=":", alpha=0.7)
    ax2.set_title("CH2 Detrending Filter (Isolates Acoustic Wavepacket)", fontsize=11, fontweight='bold')
    ax2.set_ylabel("Voltage (V)")
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper right")

    ax3.plot(t_us, np.abs(ch2_detrended), label="|CH2 Detrended Amplitude|", color="#d62728", linewidth=1.2)
    if trig_found:
        ax3.axvline(trig_t, color="blue", linestyle="--", linewidth=1.5, label=f"Trigger Start (t = {trig_t:.2f} µs)")
    ax3.axvline(peak_t, color="green", linestyle="-.", linewidth=1.5, label=f"Detected Echo Peak (t = {peak_t:.2f} µs)")
    ax3.set_title(f"Algorithm Peak Detection | Measured ToF = {calculated_tof:.3f} µs", fontsize=11, fontweight='bold')
    ax3.set_xlabel("Time (microseconds)")
    ax3.set_ylabel("Amplitude (V)")
    ax3.grid(True, linestyle="--", alpha=0.5)
    ax3.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(plot_path, dpi=150)
    plt.close()

    print(f"[SUCCESS] Scope-first debug plot saved to: file:///{plot_path.replace('\\', '/')}")

if __name__ == "__main__":
    test_correct_sync_order()
