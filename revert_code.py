# -*- coding: utf-8 -*-
import json

code = ''
with open(r'C:\Users\wille\.gemini\antigravity-ide\brain\fbea41a6-0630-415e-9764-68f6ee714fb8\.system_generated\logs\transcript.jsonl', encoding='utf-8') as f:
    for line in f:
        step = json.loads(line)
        if step.get('step_index') == 48:
            for tool_call in step.get('tool_calls', []):
                if tool_call.get('name') == 'write_to_file':
                    code = tool_call['args']['CodeContent']
                    break

old_pio = '''@rp2.asm_pio(set_init=rp2.PIO.OUT_LOW, sideset_init=rp2.PIO.OUT_LOW)
def pulse_250ns_with_scope():
    pull(block)              .side(0)       # Wait for GO signal, scope LOW
    nop()            [7]     .side(1)       # Scope HIGH, 8 cycles = 2.0 µs pre-trigger
    set(pins, 1)             .side(1)       # TX=1 (250ns pulse start), scope HIGH
    set(pins, 0)             .side(1)       # TX=0 (pulse end), scope HIGH
    set(x, 11)               .side(1)       # x=11 for 12-iteration hold loop
    label("hold_scope")
    nop()            [30]    .side(1)       # 31+1 = 32 cycles per iteration
    jmp(x_dec, "hold_scope") .side(1)       # 12 * 32 = 384 cycles = 96 µs
    nop()            [4]     .side(1)       # 5 more cycles fine-tune
    set(pins, 0)             .side(0)       # TX=0, scope LOW — all done'''

new_pio = '''@rp2.asm_pio(set_init=rp2.PIO.OUT_LOW, sideset_init=rp2.PIO.OUT_LOW)
def pulse_250ns_with_scope():
    # NOTE: With 1 sideset pin, max delay = 15 (4-bit delay field).
    # Original used [30] with no sideset (5-bit delay).
    # Restructured: 24 iterations × 16 cycles = 384 (same as 12 × 32).
    pull(block)              .side(0)       # Wait for GO signal, scope LOW
    nop()            [7]     .side(1)       # Scope HIGH, 8 cycles = 2.0 µs pre-trigger
    set(pins, 1)             .side(1)       # TX=1 (250ns pulse start), scope HIGH
    set(pins, 0)             .side(1)       # TX=0 (pulse end), scope HIGH
    set(x, 23)               .side(1)       # x=23 for 24-iteration hold loop
    label("hold_scope")
    nop()            [14]    .side(1)       # 14+1 = 15 cycles
    jmp(x_dec, "hold_scope") .side(1)       # +1 = 16 cycles per iteration, 24 × 16 = 384
    nop()            [4]     .side(1)       # 5 more cycles fine-tune
    set(pins, 0)             .side(0)       # TX=0, scope LOW — all done'''

code = code.replace(old_pio, new_pio)

old_main = '''        # Check if PC sent a serial command (e.g. "TRIG 1 2")
        events = poller.poll(0)  # 0ms non-blocking check
        if events:
            line = sys.stdin.readline().strip()'''

new_main = '''        # Check if PC sent a serial command (e.g. "TRIG 1 2")
        events = poller.poll(0)  # 0ms non-blocking check
        if events:
            char = sys.stdin.read(1)
            if char == '\\n' or char == '\\r':
                line = line_buffer.strip()
                line_buffer = ""
            else:
                line_buffer += char
                continue'''

code = code.replace(old_main, new_main)
code = code.replace('    while True:', '    line_buffer = ""\n    while True:')

with open(r'c:\Users\wille\OneDrive\Desktop\FVMFramework\examples\acoustic_imaging_with_raytracing\pico_with_ADC9226_data_acquisition\pico_scripts\pico_adc_tomography.py', 'w', encoding='utf-8') as f:
    f.write(code)

print('Success.')
