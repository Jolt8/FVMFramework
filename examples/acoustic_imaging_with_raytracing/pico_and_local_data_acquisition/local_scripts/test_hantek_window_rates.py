import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from hantek_interface_testing import Hantek6022BE

def probe_samplerate_windows():
    print("=== Testing Hantek Samplerate Divider Codes for 200µs Window ===")
    scope = Hantek6022BE()
    if not scope.connect():
        print("[X] Scope connection failed")
        return

    # Probe 0xE2 samplerate divider values:
    # 0xE2 wValue codes:
    # 0x0001: 48 MS/s -> 5.5 µs
    # 0x0002: 30 MS/s -> 8.8 µs
    # 0x0003: 16 MS/s -> 16.5 µs
    # 0x0004: 12 MS/s -> 22.0 µs
    # 0x0005: 8 MS/s  -> 33.0 µs
    # 0x0006: 4 MS/s  -> 66.0 µs
    # 0x0007: 1 MS/s  -> 264.0 µs
    # 0x0008: 500 kS/s -> 528.0 µs
    
    code_map = {
        1: (48.0, "48 MS/s"),
        2: (30.0, "30 MS/s"),
        3: (16.0, "16 MS/s"),
        4: (12.0, "12 MS/s"),
        5: (8.0, "8 MS/s"),
        6: (4.0, "4 MS/s"),
        7: (1.0, "1 MS/s"),
        8: (0.5, "500 kS/s")
    }

    try:
        scope.dev.set_interface_altsetting(0, 0)
    except Exception:
        pass

    for code, (rate_mhz, name) in code_map.items():
        try:
            scope.dev.ctrl_transfer(0x40, 0xE3, 0x0001, 0x0001, b"\x01")
            scope.dev.ctrl_transfer(0x40, 0xE2, code, 0x0000, b"\x01")
            scope.dev.ctrl_transfer(0x40, 0xE0, 0x0001, 0x0000, b"\x01")
            
            scope.dev.clear_halt(0x86)
            data = scope.dev.read(0x86, 2048, timeout=200)
            if data:
                num_samples = len(data) // 2
                window_us = num_samples / rate_mhz
                print(f"[+] Code {code:2d} ({name:10s}): Captured {len(data):4d} bytes ({num_samples:3d} samples -> {window_us:6.2f} µs window)")
        except Exception as e:
            print(f"[*] Code {code:2d} ({name:10s}): {e}")

    scope.close()

if __name__ == "__main__":
    probe_samplerate_windows()
