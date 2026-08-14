"""
Pico 2W + AD9226 12-bit ADC Tomography Firmware
================================================

Hardware-deterministic ultrasonic Time-of-Flight acquisition using:
  - PIO SM0: 250ns TX pulse (GP10) + scope trigger (GP0) at 4 MHz
  - PIO SM1: 12-bit parallel ADC capture at 10 MSPS (GP11-GP22 data)
  - Hardware PWM: Continuous 10 MHz ADC clock on GP27 (keeps ADC pipeline stable 24/7)
  - DMA: PIO RX FIFO -> RAM buffer (zero CPU intervention during capture)

Pin Assignments:
  GP0  = Scope trigger output (PIO SM0 sideset)
  GP2-5 = TX MUX channel select (4-bit)
  GP6-9 = RX MUX channel select (4-bit)
  GP10 = TX trigger output (PIO SM0 set, 250ns pulse)
  GP11 = ADC D0  (LSB, PIO SM1 in_base)
  ...
  GP22 = ADC D11 (MSB)
  GP27 = ADC CLK output (Hardware PWM, 10 MHz)
"""

import machine
import rp2
import time
import sys
import select
import array
import struct

# =============================================================================
# PIN DEFINITIONS
# =============================================================================
SCOPE_TRIGGER = 0    # GP0  — Scope trigger (PIO SM0 sideset)
TX_MUX_BASE   = 2    # GP2-5 — TX channel select
RX_MUX_BASE   = 6    # GP6-9 — RX channel select
TX_TRIGGER    = 10   # GP10 — TX pulse (PIO SM0 set)
ADC_D0_BASE   = 11   # GP11 — ADC D0 (LSB), PIO reads GP11-GP22
ADC_CLK       = 27   # GP27 — ADC clock output (Hardware PWM)
ADC_OTR       = 28   # GP28 — ADC out-of-range indicator

# =============================================================================
# DEFAULT CAPTURE PARAMETERS
# =============================================================================
DEFAULT_N_SAMPLES    = 2000       # 200 µs capture window at 10 MSPS
DEFAULT_SAMPLE_RATE  = 10_000_000 # 10 MSPS
DEFAULT_PIO_FREQ     = 125_000_000 # Run SM1 at max speed (125 MHz) to accurately wait for PWM edges

# =============================================================================
# PIO PROGRAMS
# =============================================================================

@rp2.asm_pio(set_init=rp2.PIO.OUT_LOW, sideset_init=rp2.PIO.OUT_LOW)
def pulse_250ns_with_scope():
    # 24 iterations × 16 cycles = 384 cycles.
    pull(block)              .side(0)       # Wait for GO signal, scope LOW
    nop()            [7]     .side(1)       # Scope HIGH, 8 cycles = 2.0 µs pre-trigger
    set(pins, 1)             .side(1)       # TX=1 (250ns pulse start), scope HIGH
    set(pins, 0)             .side(1)       # TX=0 (pulse end), scope HIGH
    set(x, 23)               .side(1)       # x=23 for 24-iteration hold loop
    label("hold_scope")
    nop()            [14]    .side(1)       # 14+1 = 15 cycles
    jmp(x_dec, "hold_scope") .side(1)       # +1 = 16 cycles per iteration, 24 × 16 = 384
    nop()            [4]     .side(1)       # 5 more cycles fine-tune
    set(pins, 0)             .side(0)       # TX=0, scope LOW — all done

@rp2.asm_pio(
    in_shiftdir=rp2.PIO.SHIFT_LEFT,
    autopush=True,
    push_thresh=12
)
def adc_capture_12bit():
    # PIO runs at 125 MHz (8 ns). It passively watches the 10 MHz PWM clock on GP27.
    # When GP27 goes HIGH, the AD9226 starts converting.
    # We wait for GP27 to go LOW, at which point the AD9226 data is guaranteed stable.
    wait(1, gpio, 27)               # Wait for CLK HIGH (rising edge)
    wait(0, gpio, 27)               # Wait for CLK LOW  (falling edge)
    in_(pins, 12)                   # Read 12 data pins (GP11-GP22) instantly

# =============================================================================
# HARDWARE INITIALIZATION
# =============================================================================
# --- TX Pulse PIO (SM0) ---
tx_pin = machine.Pin(TX_TRIGGER, machine.Pin.OUT)
scope_pin = machine.Pin(SCOPE_TRIGGER, machine.Pin.OUT)
sm0 = rp2.StateMachine(
    0,                          
    pulse_250ns_with_scope,
    freq=4_000_000,             
    set_base=tx_pin,            
    sideset_base=scope_pin      
)
sm0.active(1)

# --- ADC Clock (Hardware PWM on GP27) ---
# Running this continuously ensures the pipelined AD9226 NEVER freezes its internal biases.
adc_clk_pwm = machine.PWM(machine.Pin(ADC_CLK))
adc_clk_pwm.freq(DEFAULT_SAMPLE_RATE)
adc_clk_pwm.duty_u16(32768)  # 50% duty cycle

# --- ADC Capture PIO (SM1) ---
adc_d0_pin = machine.Pin(ADC_D0_BASE, machine.Pin.IN)
for i in range(12):
    machine.Pin(ADC_D0_BASE + i, machine.Pin.IN)

sm1 = rp2.StateMachine(
    1,                          
    adc_capture_12bit,
    freq=DEFAULT_PIO_FREQ,      
    in_base=adc_d0_pin
)

# --- MUX GPIO ---
tx_mux_pins = [machine.Pin(TX_MUX_BASE + i, machine.Pin.OUT) for i in range(4)]
rx_mux_pins = [machine.Pin(RX_MUX_BASE + i, machine.Pin.OUT) for i in range(4)]

# =============================================================================
# DMA CONFIGURATION
# =============================================================================
# DREQ for PIO0 SM1 RX is ALWAYS 3 on both RP2040 and RP2350 (Datasheet 2.5.3.1)
DREQ_PIO0_RX1 = 3

dma = rp2.DMA()

n_samples = DEFAULT_N_SAMPLES
sample_rate_hz = DEFAULT_SAMPLE_RATE
capture_buf = array.array('I', [0] * n_samples)
send_buf = bytearray(n_samples * 2)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def set_mux(pins, channel):
    for i in range(4):
        pins[i].value((channel >> i) & 1)

def fire_tx_pulse():
    sm0.put(1)

def reallocate_buffers():
    global capture_buf, send_buf
    capture_buf = array.array('I', [0] * n_samples)
    send_buf = bytearray(n_samples * 2)

def reconfigure_adc(new_rate_hz, new_n_samples):
    global n_samples, sample_rate_hz
    n_samples = new_n_samples
    sample_rate_hz = new_rate_hz
    adc_clk_pwm.freq(new_rate_hz)
    reallocate_buffers()

def capture_adc():
    # Drain stale FIFO data
    sm1.active(0)
    while sm1.rx_fifo() > 0:
        sm1.get()
    
    # Configure DMA transfer
    dma.config(
        read=sm1,
        write=capture_buf,
        count=n_samples,
        ctrl=dma.pack_ctrl(
            size=2,
            treq_sel=DREQ_PIO0_RX1,
            inc_read=False,
            inc_write=True
        ),
        trigger=False
    )
    
    # Start DMA and PIO
    dma.active(1)
    sm1.active(1)
    
    # Fire TX pulse
    fire_tx_pulse()
    
    # Wait for completion (200µs @ 10 MSPS)
    timeout_us = max(n_samples * 2, 1000)
    t_start = time.ticks_us()
    while dma.active():
        if time.ticks_diff(time.ticks_us(), t_start) > timeout_us:
            break
            
    sm1.active(0)
    return capture_buf

def pack_and_send(buf, n):
    for i in range(n):
        val = buf[i] & 0x0FFF
        send_buf[i*2]     = val & 0xFF
        send_buf[i*2 + 1] = (val >> 8) & 0xFF
    sys.stdout.write(f"DATA {n} {sample_rate_hz}\n")
    sys.stdout.buffer.write(send_buf[:n*2])
    sys.stdout.write("END\n")

# =============================================================================
# MAIN LOOP
# =============================================================================
def main():
    transducer_send_receive_ordering = [(1, 2)]
    
    print("Pico 2W + AD9226 ADC Tomography Controller Ready (PWM CLOCK v2).")
    
    poller = select.poll()
    poller.register(sys.stdin, select.POLLIN)
    line_buffer = ""
    
    while True:
        command_ready = None
        while poller.poll(0):
            char = sys.stdin.read(1)
            if char == '\n' or char == '\r':
                if line_buffer:
                    command_ready = line_buffer.strip()
                    line_buffer = ""
                    break
            else:
                line_buffer += char
                
        if command_ready:
            line = command_ready
            
            if line == "PING":
                print("PONG")
            
            elif line.startswith("TRIG"):
                parts = line.split()
                if len(parts) >= 3:
                    set_mux(tx_mux_pins, int(parts[1]))
                    set_mux(rx_mux_pins, int(parts[2]))
                    time.sleep_us(10)
                    
                    buf = capture_adc()
                    pack_and_send(buf, n_samples)
                else:
                    print("ERR bad TRIG format")
            
            elif line.startswith("CONF"):
                parts = line.split()
                if len(parts) >= 3:
                    new_rate = int(parts[1])
                    new_n = int(parts[2])
                    if 100_000 <= new_rate <= 65_000_000 and 10 <= new_n <= 50000:
                        reconfigure_adc(new_rate, new_n)
                        print(f"CONF_ACK {sample_rate_hz} {n_samples}")
                else:
                    print("ERR bad CONF format")
            
            elif line == "STATUS":
                print(f"STATUS rate={sample_rate_hz} n={n_samples} otr=0")
            
            else:
                print(f"ERR unknown command: {line}")
            
            continue
        
        # --- STANDALONE CONTINUOUS-FIRE MODE ---
        for tx, rx in transducer_send_receive_ordering:
            set_mux(tx_mux_pins, tx)
            set_mux(rx_mux_pins, rx)
            fire_tx_pulse()
            time.sleep_us(250)

if __name__ == "__main__":
    main()
