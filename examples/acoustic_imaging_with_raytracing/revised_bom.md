# Revised BOM — 16-Element 2 MHz Ultrasonic Transducer Array

## System Architecture

```
                          ┌─────────────┐
                          │  Pico 2     │
                          │ (Controller)│
                          └──┬───┬───┬──┘
                    4 GPIO   │   │   │  4 GPIO
                (TX addr)    │   │   │  (RX addr)
              ┌──────────────┘   │   └───────────────┐
              ▼                  │                    ▼
     ┌────────────────┐         │          ┌────────────────┐
     │ HC4067 #1      │         │          │ HT6022BE Scope │
     │ (TX gate mux)  │         │          │ (digitizer)    │
     │ COM pin → 5V   │         │          └───────▲────────┘
     └───┬──┬──┬──────┘         │                  │
         │  │  │ (16 gate lines)│          ┌───────┴────────┐
         ▼  ▼  ▼                │          │ AD828 board    │
     ┌──────────────┐           │          │ (RX amplifier) │
     │ 16× 2N7000   │           │          └───────▲────────┘
     │ (TX switches)│           │                  │
     └───┬──┬──┬────┘           │          ┌───────┴────────┐
         │  │  │                │          │ Bandpass Filter │
         │  │  │  common TX bus │          │ (1–3 MHz)      │
     ┌───┴──┴──┴────┐          │          └───────▲────────┘
     │   TC4420      │          │                  │
     │  (TX pulser)  │          │          ┌───────┴────────┐
     └───────────────┘          │          │ HC4067 #2      │
                                │          │ (RX mux)       │
                                │          └───▲──▲──▲──────┘
                                │              │  │  │
         ┌──────────────────────┘              │  │  │
         │           ┌────────────────────────┘  │  │
         │           │              ┌────────────┘  │
     ┌───┴───────────┴──────────────┴───────────────┴───┐
     │  16× 1N4148 pairs (T/R protection per element)   │
     └───┬───────────┬──────────────┬───────────────┬───┘
         │           │              │               │
     ┌───┴───┐   ┌───┴───┐     ┌───┴───┐       ┌───┴───┐
     │ XDCR  │   │ XDCR  │     │ XDCR  │  ...  │ XDCR  │
     │  #1   │   │  #2   │     │  #3   │       │  #16  │
     └───────┘   └───────┘     └───────┘       └───────┘
```

> [!IMPORTANT]
> **Key design principle:** TX and RX paths are separated. The TX HC4067 only switches MOSFET gates (logic-level, zero current). Discrete MOSFETs handle TX pulse switching (high current). The RX HC4067 routes small analog echo signals. Diode limiters protect the RX path during transmit.
>
> **Independent TX/RX selection:** Each HC4067 has its own 4 address lines (8 GPIO total), allowing you to transmit on one element and receive on any other — enabling both pulse-echo and through-transmission modes.

---

## BOM by Subsystem

### 1. Transducers

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 1 | 2 MHz ultrasonic transducer | 16 | Module/element | Using 16 of the 20 purchased. Ensure you know the impedance (~50Ω typical) for matching. Keep the remaining 4 as spares. |

---

### 2. Controller

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 2 | Raspberry Pi Pico 2 | 1 | **Dev board** | Buy the board with headers pre-soldered if you want easy breadboarding. Provides 3.3V logic and 30 GPIO pins. Uses 8 GPIO for mux addressing + 1 for TC4420 trigger + 1 for scope trigger = 10 total, leaving 20 free. |

> [!TIP]
> The Pico 2 uses 3.3V logic. The HC4067 and 2N7000 gates work fine at 3.3V (2N7000 threshold is ~2.1V). For the TC4420 input (4.5V logic threshold), you'll need a **level shifter** on the trigger line — see the level shifter entry below.

---

### 3. TX Pulse Generation

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 4 | TC4420 | 1 | **Bare DIP-8 IC** | Single MOSFET driver. Generates the TX pulse on a common bus. Supply voltage sets your pulse amplitude (up to 18V). Start at 5–12V and increase if SNR is too low. |
| 5 | 10Ω gate resistor | 1 | Through-hole | Series resistor on TC4420 output to damp ringing. |
| 6 | 0.1µF ceramic capacitor | 1 | Through-hole | Bypass cap close to TC4420 VDD pin. |
| 7 | 10µF electrolytic capacitor | 1 | Through-hole | Bulk decoupling for TC4420 supply. |

> [!NOTE]
> **On TX voltage:** The TC4420 at 12V will produce a short unipolar pulse. This is fine for short-range experiments in water or thin-material NDT. If you find the echo signal is too weak, you can upgrade later to a MOSFET half-bridge with a 50–100V supply, but **start simple**.

---

### 4. TX Element Selection (HC4067 #1 + MOSFETs)

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 3 | CD74HC4067 (TX gate mux) | 1 | **Breakout board** | Used as a 1-to-16 demux. Common pin tied to 5V; selected output goes HIGH to open one MOSFET gate. Same breakout boards as the RX mux. ~$1–2. |
| 8 | 2N7000 N-channel MOSFET | 16 | **Bare TO-92** | One per transducer. Acts as a switch between the TX bus and each element. V_DS = 60V, I_D = 200mA, R_DS(on) ~5Ω. TO-92 package is breadboard/perfboard friendly. |
| 9 | 10kΩ resistor | 16 | Through-hole | Gate pulldown for each 2N7000 — ensures OFF state when not driven by the mux. |
| 10 | 100Ω resistor | 16 | Through-hole | Series between TX bus and transducer (through the MOSFET drain) to limit ringing and current. Adjust value based on your transducer impedance. |

**How it works:**
- The Pico sets 4 address lines on HC4067 #1 → selected output goes to 5V → turns on one 2N7000 gate.
- The TC4420 fires a pulse → current flows through the selected MOSFET → into the selected transducer.
- All other MOSFET gates are pulled LOW by their 10kΩ resistors (the HC4067's non-selected outputs are high-impedance, so the pulldowns hold them OFF).
- The RX HC4067 (#2) can independently select a different element for receiving.

---

### 5. T/R Protection

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 11 | 1N4148 fast diode | 32 | **Bare through-hole** | **2 per element** in a back-to-back clamp configuration (one to VCC, one to GND) at the RX mux input. Clamps the TX pulse to ~±0.7V to protect the RX amplifier. ~$0.02 each. |

**Wiring per element (protection clamp):**
```
          Transducer element
                  │
                  ├───────── to TX MOSFET drain
                  │
       ┌──────────┤ (RX tap)
       │          │
    ┌──┴──┐    ┌──┴──┐
    │D1   │    │D2   │
    │1N4148│    │1N4148│
    └──┬──┘    └──┬──┘
       │          │
      VCC        GND
       │
       └──── to RX mux input
```

> [!NOTE]
> Alternatively, you can buy **20× BAV99** (dual series diode in SOT-23, ~$0.05 each) which contain two diodes per package and have lower parasitic capacitance (~1.5 pF vs ~4 pF). However, SOT-23 is harder to breadboard — only worth it if you're doing a PCB.

---

### 6. RX Multiplexing (HC4067 #2)

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 12 | CD74HC4067 (RX mux) | 1 | **Breakout board** | Same type of board as the TX gate mux. Routes the selected transducer's echo signal to the amplifier chain. |

> [!TIP]
> Both HC4067 boards are identical hardware — the only difference is how they're wired (TX mux: common=5V, outputs to MOSFET gates; RX mux: common=amplifier input, outputs from transducer RX taps). If you later move to a PCB, upgrade to **ADG1606** (TSSOP-24, ~$6, 4Ω on-resistance) for better RX signal integrity.

**Pico GPIO usage for RX mux:** 4 address lines (S0–S3) = **4 GPIO pins** (separate from the TX mux's 4 lines).

**Total GPIO for muxes:** 8 pins (4 TX + 4 RX) = independent TX/RX element selection.

---

### 7. RX Amplification

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 13 | AD828 dual op-amp | 1 | **Pre-built board** | Your board already has gain-setting resistors on it — no need to buy separate ones. Check the board's documentation for the preset gain value. If the gain is wrong for your application, you can swap the on-board resistors later. |
| 15 | 0.1µF ceramic capacitor | 2 | Through-hole | Bypass caps on ±V supply pins of the AD828. Ceramic is preferred over electrolytic here — better high-frequency decoupling at 2 MHz. |
| 17 | 10nF DC-blocking cap | 1 | Through-hole/ceramic | AC coupling at the input to block DC offset from the RX mux. If your board already has input coupling caps, skip this. |

> [!IMPORTANT]
> The AD828 needs a **dual supply** (e.g., ±5V) for best performance, or a single supply (10V) with a mid-rail bias. Check whether your board handles this or if you need the ICL7660 voltage inverter (see Power Supply section).

---

### 8. Bandpass Filter (1–3 MHz)

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 18 | 10µH inductor | 2 | Through-hole axial | Series elements of a simple LC bandpass. |
| 19 | 1nF ceramic capacitor (C0G/NP0) | 2 | Through-hole | Shunt elements. These values give a center frequency near 2 MHz. Use **C0G/NP0** dielectric for stable, low-loss RF behavior — do NOT use X7R or Y5V. |

**Simple LC bandpass topology:**
```
Input ──┤├── L1 ──┬── L2 ──┤├── Output
       10nF       │       10nF
              ┌───┴───┐
              │  1nF  │
              └───┬───┘
                  │
                 GND
```

> [!NOTE]
> This is a rough starting point. You will likely need to tune L and C values on the bench using your scope to center the passband on your transducer's actual resonant frequency. Buy an **assortment** of NP0 caps (100pF–10nF) and inductors (1µH–100µH) to have tuning flexibility.

---

### 9. Digitizer

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 20 | HT6022BE 20 MHz USB scope | 1 | **Complete device** | 48 MSa/s, 8-bit. Connect to a PC; trigger from a Pico GPIO output. Software: OpenHantek or Hantek official. Essential for debugging waveforms, checking pulse shape, measuring timing, and diagnosing issues. |

> [!NOTE]
> The scope is your most important debugging tool. For automated multi-channel scanning later, you could add a dedicated ADC, but the scope should remain in your toolkit permanently.

---

### 10. Level Shifting

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 21 | 4-channel bidirectional level shifter (3.3V ↔ 5V) | 1 | **Breakout board** | BSS138-based 4-ch module, ~$1. Needed only for the **TC4420 trigger line** (requires ≥4.5V logic HIGH). The HC4067 address lines and 2N7000 gates work fine at 3.3V directly from the Pico. |

> [!NOTE]
> The 2N7000's gate threshold is ~2.1V and the HC4067 accepts 3.3V logic when powered from 3.3V or 5V (V_IH = 0.7×VCC). So only the TC4420 trigger needs level shifting — 1 channel of the 4-ch module, leaving 3 spare channels.

---

### 11. Power Supply

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 22 | 5V 2A USB power supply | 1 | Wall adapter / USB | Powers the Pico (via USB), logic ICs, and RX amplifier (single supply mode). |
| 23 | Adjustable bench power supply (0–30V) | 1 | Bench unit | For the TX pulse voltage. Start at 5V, increase to 12–18V as needed. **You may already own one.** |
| 24 | ICL7660 or TC1044 (voltage inverter) | 1 | **Bare DIP-8 IC** | Converts +5V to −5V to create a ±5V dual supply for the AD828. Only needed if your AD828 board doesn't handle dual supply on its own. ~$1. |
| 25 | 0.1µF ceramic capacitor | 10 | Through-hole | General-purpose decoupling. Place one near VCC of every IC. Ceramic is fine for all decoupling in this system — no electrolytics needed. |

---

### 12. Cabling & Connectors

| # | Component | Qty | Form Factor | Notes |
|---|-----------|-----|-------------|-------|
| 27 | RG174 coax cable | ~10 m | Bulk cable | Thin 50Ω coax suitable for 2 MHz. Cut to length for each transducer. RG174 is flexible and easy to work with. |
| 28 | SMA connectors (panel mount + cable) | 20 pairs | Connectors | If you want proper detachable transducer connections. Otherwise, just solder directly for prototyping. **Optional for Phase 1.** |
| 29 | Breadboard (full size) | 2 | Board | For prototyping. One for TX, one for RX. |
| 30 | Jumper wire kit | 1 | Assortment | Male-male, male-female, female-female. |
| 31 | Perfboard (9×15 cm) | 2 | Board | For when you move off breadboard. **Phase 2.** |

---

## Complete Summary Table

| # | Component | Qty | IC / Board | Est. Unit Cost | Est. Total |
|---|-----------|-----|------------|---------------|------------|
| 1 | 2 MHz ultrasonic transducer | 16 | Module | varies | varies |
| 2 | Raspberry Pi Pico 2 | 1 | Board | $5 | $5 |
| 3 | CD74HC4067 (TX gate mux) | 1 | Breakout board | $1.50 | $1.50 |
| 4 | TC4420 MOSFET driver | 1 | Bare DIP IC | $2.50 | $2.50 |
| 5 | 2N7000 N-ch MOSFET | 16 | Bare TO-92 | $0.15 | $2.40 |
| 6 | CD74HC4067 (RX mux) | 1 | Breakout board | $1.50 | $1.50 |
| 7 | AD828 dual op-amp | 1 | Pre-built board | $4 | $4 |
| 8 | 1N4148 fast diode | 32 | Bare through-hole | $0.02 | $0.64 |
| 9 | 4-ch level shifter | 1 | Breakout board | $1 | $1 |
| 10 | ICL7660 voltage inverter | 1 | Bare DIP IC | $1 | $1 |
| 11 | 10kΩ resistor | 16 | Through-hole | — | $1 (pack) |
| 12 | 100Ω resistor | 16 | Through-hole | — | $1 (pack) |
| 13 | 10Ω resistor | 1 | Through-hole | — | — |
| 14 | 0.1µF ceramic cap | 12 | Through-hole | — | $1 (pack) |
| 15 | 10nF DC-blocking cap | 1 | Through-hole | — | — |
| 16 | 1nF NP0 ceramic cap | 2 | Through-hole | — | — |
| 17 | 10µH inductor | 2 | Through-hole axial | $0.20 | $0.40 |
| 18 | HT6022BE USB scope | 1 | Complete device | $45 | $45 |
| 19 | RG174 coax cable | 8 m | Bulk | $0.50/m | $4 |
| 20 | Breadboard (full size) | 2 | Board | $5 | $10 |
| 21 | Jumper wire kit | 1 | Assortment | $5 | $5 |
| | | | | **Subtotal (excl. transducers, bench supply, scope):** | **~$37** |
| | | | | **Subtotal (with scope):** | **~$82** |

---

## Phased Build Plan

### Phase 1 — Single-Element Proof of Concept
Build and verify with **1 transducer** (no mux needed):
- Pico → TC4420 → transducer → 1N4148 clamp → AD828 board → HT6022BE
- **Goal**: See a TX pulse and a received echo on the scope.

### Phase 2 — Add Both Muxes
- Add 2× HC4067 breakout boards (one for TX gate selection, one for RX).
- Wire up 16× 2N7000 MOSFETs with gate pulldowns.
- Wire up all 16 T/R protection diode pairs.
- Connect 8 GPIO lines from Pico (4 per mux).
- **Goal**: Independently select TX and RX elements. Verify pulse-echo on any element AND through-transmission between any pair.

### Phase 3 — Full Array Commissioning
- Add bandpass filter, verify signal quality across all 16 channels.
- Write Pico firmware for automated scanning sequences.
- Characterize SNR, crosstalk, and timing for each element pair.

---

## What Changed from Your Original List

| Original Component | Change | Reason |
|---|---|---|
| 1× CD74HC4067 | Now using **2× HC4067 breakout boards** | One as TX gate demux, one as RX signal mux. Independent TX/RX element selection via 8 GPIO lines. |
| Multiple TC4420s | Reduced to 1 | Only 1 needed since TX is multiplexed through MOSFETs. |
| AD828 as bare IC | Using your **pre-built AD828 board** | No separate gain resistors needed — they're on the board. |
| Electrolytic caps | **Replaced with ceramics throughout** | Better high-frequency performance at 2 MHz, and you already have a ceramic cap kit. |
| HT6022BE | **Kept as primary digitizer** | Best debugging tool; AD9220 upgrade deferred to future if needed. |

## What Was Added

| Added Component | Reason |
|---|---|
| 2N7000 MOSFETs | TX element switches — handle pulse current that the HC4067 cannot. |
| Level shifter (1× 4-ch) | TC4420 trigger line needs 5V logic; only 1 channel used. |
| ICL7660 voltage inverter | Create ±5V supply for the AD828 (if your board needs it). |
| Bandpass filter (L + C) | Reject out-of-band noise before amplification. |
| Decoupling capacitors (ceramic) | Essential for stable IC operation — not optional. |
| DC-blocking capacitor | Prevent DC offset from saturating amplifier input. |
