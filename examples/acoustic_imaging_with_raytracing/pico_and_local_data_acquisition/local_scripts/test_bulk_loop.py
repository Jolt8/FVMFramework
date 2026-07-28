import sys
import time
import libusb_package
import usb.core
import usb.util

def test_hantek_continuous():
    print("=== Testing Hantek 6022BE Continuous Streaming Mode ===")
    backend = libusb_package.get_libusb1_backend()
    dev = usb.core.find(idVendor=0x04B5, idProduct=0x6022, backend=backend)

    if not dev:
        print("[X] Device not found!")
        return

    try:
        dev.set_configuration()
        usb.util.claim_interface(dev, 0)
        print("[+] Claimed USB Interface 0")
    except Exception as e:
        print(f"[X] Claim interface error: {e}")
        return

    try:
        dev.set_interface_altsetting(0, 0)
    except Exception:
        pass

    # OpenHantek 6022BE startup control sequence:
    # 0xE3 (Gain): wValue=0x0001, wIndex=0x0001
    # 0xE2 (Samplerate): wValue=0x0001
    # 0xE4 (Trigger mode): wValue=0x0001 (Continuous)
    # 0xE0 (Start): wValue=0x0001
    dev.ctrl_transfer(0x40, 0xE3, 0x0001, 0x0001, b"")
    dev.ctrl_transfer(0x40, 0xE2, 0x0001, 0x0000, b"")
    dev.ctrl_transfer(0x40, 0xE4, 0x0001, 0x0000, b"")
    dev.ctrl_transfer(0x40, 0xE0, 0x0001, 0x0000, b"")
    time.sleep(0.05)

    dev.clear_halt(0x86)

    raw_bytes = bytearray()
    t0 = time.time()
    
    # Read bulk endpoint in loop
    for i in range(100):
        try:
            chunk = dev.read(0x86, 2048, timeout=50)
            if chunk:
                raw_bytes.extend(chunk)
        except Exception:
            pass

    t_elapsed = time.time() - t0
    num_samples = len(raw_bytes) // 2
    duration_us = num_samples / 16.0

    print(f"[+] Read {len(raw_bytes)} bytes in {t_elapsed*1000.0:.1f} ms")
    print(f"[+] Total window duration captured: {duration_us:.2f} µs!")

    try:
        usb.util.release_interface(dev, 0)
    except Exception:
        pass

if __name__ == "__main__":
    test_hantek_continuous()
