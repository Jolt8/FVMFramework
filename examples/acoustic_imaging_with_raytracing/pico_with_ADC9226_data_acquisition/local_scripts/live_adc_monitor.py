"""
Live ADC Monitor
================
Continuously triggers the ADC and prints the average raw ADC code and voltage.
Use this to debug the AD9226 hardware.
"""
import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pico_adc_interface import PicoADCInterface

def main():
    print("=== LIVE ADC MONITOR ===")
    adc = PicoADCInterface()
    if not adc.connect():
        return
        
    print("Press Ctrl+C to stop.\n")
    try:
        while True:
            volts, rate, dt = adc.trigger_and_capture(tx=1, rx=2)
            if volts is not None:
                raw_codes = (volts / adc.vref * 2048.0) + 2048.0
                mean_code = np.mean(raw_codes)
                mean_volt = np.mean(volts)
                std_code = np.std(raw_codes)
                
                print(f"Mean Code: {mean_code:4.0f} | Mean Volt: {mean_volt:+.3f} V | Noise (Std): {std_code:.1f} LSBs")
            time.sleep(0.1)
    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        adc.close()

if __name__ == "__main__":
    main()
