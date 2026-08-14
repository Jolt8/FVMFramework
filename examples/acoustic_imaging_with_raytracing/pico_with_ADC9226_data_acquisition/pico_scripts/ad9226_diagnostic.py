"""
AD9226 Standalone Test v2 — Non-blocking with verbose progress
"""
import machine
import rp2
import time

ADC_D0_BASE = 11
ADC_CLK     = 1

print("=" * 60)
print("  AD9226 STANDALONE TEST v2")
print("=" * 60)

# --- Step 1: Start clock ---
print("\n[1] Starting 10 MHz clock on GP1...")
clk = machine.PWM(machine.Pin(ADC_CLK))
clk.freq(10_000_000)
clk.duty_u16(32768)
print(f"    OK. Reported freq={clk.freq()} Hz")

# --- Step 2: Wait for ADC pipeline ---
print("[2] Waiting 100ms for ADC pipeline to fill...")
time.sleep_ms(100)
print("    OK.")

# --- Step 3: Read data pins directly (no PIO) ---
print("[3] Direct GPIO reads (sanity check):")
data_pins = [machine.Pin(ADC_D0_BASE + i, machine.Pin.IN) for i in range(12)]
for trial in range(3):
    word = 0
    for i in range(12):
        word |= (data_pins[i].value() << i)
    print(f"    Trial {trial}: code={word:4d} (0x{word:03X})")
    time.sleep_ms(10)

# --- Step 4: PIO capture with timeout ---
print("[4] Setting up PIO...")

@rp2.asm_pio(in_shiftdir=rp2.PIO.SHIFT_LEFT, autopush=True, push_thresh=12)
def adc_read():
    wait(1, gpio, 1)
    wait(0, gpio, 1)
    in_(pins, 12)

sm = rp2.StateMachine(0, adc_read, freq=125_000_000,
                       in_base=machine.Pin(ADC_D0_BASE, machine.Pin.IN))
print("    PIO configured.")

print("[5] Activating PIO and checking FIFO (2 second timeout)...")
sm.active(1)

# Non-blocking: poll FIFO with timeout
deadline = time.ticks_ms() + 2000
fifo_count = 0
while time.ticks_diff(time.ticks_ms(), deadline) < 0:
    fifo_count = sm.rx_fifo()
    if fifo_count > 0:
        break
    time.sleep_ms(10)

if fifo_count == 0:
    print("    ** FIFO EMPTY after 2s! PIO is stuck on wait(). **")
    print("    ** GP1 clock is not reaching the PIO. **")
    sm.active(0)
    
    # Fallback: try reading WITHOUT clock sync
    print("\n[6] Fallback: PIO reads WITHOUT clock sync (free-running)...")
    
    @rp2.asm_pio(in_shiftdir=rp2.PIO.SHIFT_LEFT, autopush=True, push_thresh=12)
    def adc_freerun():
        in_(pins, 12)
        nop() [31]
        nop() [31]
    
    sm2 = rp2.StateMachine(1, adc_freerun, freq=125_000_000,
                            in_base=machine.Pin(ADC_D0_BASE, machine.Pin.IN))
    sm2.active(1)
    time.sleep_ms(50)
    
    print("    Free-running PIO samples:")
    for i in range(20):
        if sm2.rx_fifo() > 0:
            raw = sm2.get()
            code = raw & 0x0FFF
            volts = (code - 2048.0) / 2048.0
            print(f"      Sample {i:2d}: code={code:4d} | volts={volts:+.4f} V")
        else:
            print(f"      Sample {i:2d}: FIFO empty")
        time.sleep_ms(10)
    sm2.active(0)

else:
    print(f"    FIFO has {fifo_count} entries! PIO is capturing data!")
    
    # Read 20 samples
    print("\n[6] Reading 20 samples:")
    for i in range(20):
        if sm.rx_fifo() > 0:
            raw = sm.get()
            code = raw & 0x0FFF
            volts = (code - 2048.0) / 2048.0
            print(f"    Sample {i:2d}: code={code:4d} (0x{code:03X}) | volts={volts:+.4f} V")
        else:
            # Wait briefly for next sample
            time.sleep_us(10)
            if sm.rx_fifo() > 0:
                raw = sm.get()
                code = raw & 0x0FFF
                volts = (code - 2048.0) / 2048.0
                print(f"    Sample {i:2d}: code={code:4d} (0x{code:03X}) | volts={volts:+.4f} V")
            else:
                print(f"    Sample {i:2d}: FIFO empty")
    
    sm.active(0)

# --- Step 7: Summary ---
print("\n[7] Reading all 12 pins one more time (clock still running):")
word = 0
for i in range(12):
    word |= (data_pins[i].value() << i)
print(f"    Final code: {word:4d} (0x{word:03X})")

clk.deinit()
print("\n" + "=" * 60)
print("  TEST COMPLETE")
print("=" * 60)
