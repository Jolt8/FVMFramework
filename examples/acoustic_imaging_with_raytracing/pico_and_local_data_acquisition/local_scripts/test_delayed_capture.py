import sys
import time
import os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from tomography_master_acquisition import TomographyAcquisitionMaster

def test_full_tof_calculation(pre_delay_us=50.0):
    print(f"=== Sandboxing Precision ToF Calculation with {pre_delay_us} µs Pre-Delay ===")
    master = TomographyAcquisitionMaster()
    if not master.connect():
        print("[X] Connection failed")
        return

    # Trigger Pico
    master.pico.trigger_pair(1, 2)
    time.sleep(pre_delay_us / 1e6)

    ch1, ch2 = master.scope.capture_waveform(num_samples=2048)
    if ch1 is not None and ch2 is not None:
        dt_us = 1.0 / 16.0
        t_us = np.arange(len(ch2)) * dt_us + pre_delay_us
        
        # Locate acoustic pulse peak on CH2
        v_rx = ch2 - np.mean(ch2)
        peak_sample_idx = np.argmax(np.abs(v_rx))
        tof_us = pre_delay_us + (peak_sample_idx * dt_us)
        
        print(f"[SUCCESS!] CH2 Acoustic Wavepacket Peak detected at sample #{peak_sample_idx}!")
        print(f"[SUCCESS!] Calculated True Time-of-Flight (ToF): {tof_us:.3f} µs (Target ~68 µs)")

    master.scope.close()
    master.pico.close()

if __name__ == "__main__":
    test_full_tof_calculation(pre_delay_us=50.0)
