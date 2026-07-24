import machine
import rp2
import time
import sys
import select

# --- PIN DEFINITIONS ---
TX_MUX_BASE = 2     # GP2, GP3, GP4, GP5 (4-bit TX Channel Selection)
RX_MUX_BASE = 6     # GP6, GP7, GP8, GP9 (4-bit RX Channel Selection)
TX_TRIGGER = 10     # GP10 (Driven by PIO state machine for 250ns pulse)
SCOPE_TRIGGER = 11  # GP11 (Scope trigger output)

# --- PIO PROGRAM FOR EXACT 250ns HARDWARE PULSE ---
# The PIO state machine runs independently in hardware at 4 MHz.
# 1 clock cycle at 4 MHz = 1 / 4,000,000 sec = exactly 250 nanoseconds.
@rp2.asm_pio(set_init=rp2.PIO.OUT_LOW)
def pulse_250ns():
    pull(block)     # Wait until MicroPython puts a value into the FIFO
    set(pins, 1)    # Set TX pin HIGH (lasts exactly 1 clock cycle = 250ns)
    set(pins, 0)    # Set TX pin LOW

# --- INITIALIZE PIO STATE MACHINE ---
tx_pin = machine.Pin(TX_TRIGGER, machine.Pin.OUT)

# Setting freq=4_000_000 configures the PIO clock divider automatically
sm = rp2.StateMachine(0, pulse_250ns, freq=4_000_000, set_base=tx_pin)
sm.active(1)

# --- INITIALIZE GENERAL GPIO PINS ---
tx_mux_pins = [machine.Pin(TX_MUX_BASE + i, machine.Pin.OUT) for i in range(4)]
rx_mux_pins = [machine.Pin(RX_MUX_BASE + i, machine.Pin.OUT) for i in range(4)]
scope_trigger_pin = machine.Pin(SCOPE_TRIGGER, machine.Pin.OUT, value=0)

def set_mux(pins, channel):
    """Sets a 4-bit multiplexer channel given a list of 4 GPIO pins."""
    for i in range(4):
        pins[i].value((channel >> i) & 1)

def send_pulse():
    """
    Triggers the PIO hardware to fire an exact 250ns pulse.
    MicroPython doesn't sleep here; it just signals the PIO hardware block.
    """
    sm.put(1)

def trigger_pair(tx_chan, rx_chan):
    """Sets multiplexers, raises scope trigger, and fires 250ns ultrasonic pulse."""
    set_mux(tx_mux_pins, tx_chan)
    set_mux(rx_mux_pins, rx_chan)
    scope_trigger_pin.value(1)
    send_pulse()
    time.sleep_us(50)
    scope_trigger_pin.value(0)

def main():
    transducer_send_receive_ordering = [
        (1, 2), (2, 1),
        (3, 4), (4, 3),
        (5, 6), (6, 5),
        (7, 8), (8, 7)
    ]

    print("Pico 2W Tomography Controller Ready (USB Serial + Standalone Mode).")

    poller = select.poll()
    poller.register(sys.stdin, select.POLLIN)

    while True:
        # Check if PC sent a serial command (e.g. "TRIG 1 2")
        events = poller.poll(10)  # 10ms timeout
        if events:
            line = sys.stdin.readline().strip()
            if line.startswith("TRIG"):
                parts = line.split()
                if len(parts) >= 3:
                    tx_c = int(parts[1])
                    rx_c = int(parts[2])
                    trigger_pair(tx_c, rx_c)
                    print(f"ACK {tx_c} {rx_c}")
                    continue

        # Standalone loop mode if no PC command
        for send_transducer, receive_transducer in transducer_send_receive_ordering:
            trigger_pair(send_transducer, receive_transducer)
            time.sleep_ms(10)

if __name__ == "__main__":
    main()
