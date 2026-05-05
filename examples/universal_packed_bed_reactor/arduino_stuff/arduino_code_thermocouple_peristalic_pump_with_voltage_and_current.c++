// Merged Arduino code: 5 thermocouples + peristaltic pump + AC voltage/current
// The pump continues stepping during ALL blocking operations (RMS sampling,
// thermocouple retries, etc.) via pumpDuringWait() helpers.
//
// Hardware summary:
//   - 5x MAX31855 thermocouples on SPI (CS pins 2-6)
//   - TB6600 stepper driver (PUL 9, DIR 7, ENA 8)
//   - ZMPT101B voltage sensor on A0
//   - SCT-013 current sensor via ADS1115 on I2C (differential A0/A1)

#include <SPI.h>
#include <Wire.h>
#include <Adafruit_MAX31855.h>
#include <Adafruit_ADS1X15.h>

// =============================================================
//  THERMOCOUPLE SETUP
// =============================================================
#define CS1 2
#define CS2 3
#define CS3 4
#define CS4 5
#define CS5 6

Adafruit_MAX31855 tc1(CS1);
Adafruit_MAX31855 tc2(CS2);
Adafruit_MAX31855 tc3(CS3);
Adafruit_MAX31855 tc4(CS4);
Adafruit_MAX31855 tc5(CS5);

Adafruit_MAX31855* thermocouples[] = {&tc1, &tc2, &tc3, &tc4, &tc5};
const int numSensors = 5;

// Non-blocking temperature read interval
unsigned long previousTempMillis = 0;
const long tempInterval = 500; // ms between temperature reads

// =============================================================
//  STEPPER MOTOR (PERISTALTIC PUMP) SETUP
// =============================================================
#define DIR_PIN 7
#define ENA_PIN 8
#define PUL_PIN 9

unsigned long previousStepMicros = 0;

double target_flow_rate_ml_per_min = 1.0; // ← change me!
double target_rpm = -0.07234923187935079 + 3.2770887868923912 * target_flow_rate_ml_per_min;

long pulses_per_revolution = 6400;
// check TB6600 for the pulse/rev, the setting resulting in the highest pulse/rev should be used 
// unless you need more head pressure or need to overcome a higher pressure drop
// I'm not sure if you would need to do another interpolation if you change the pulse/rev
// since this is a positive displacement pump where the relationship between the flow rate and rpm is very linear, 
// I don't think a new interpolation would be required
double pulses_per_second = target_rpm * pulses_per_revolution / 60.0;
unsigned long stepIntervalMicros = 1000000UL / (unsigned long)pulses_per_second;

bool motorDirection = true; // true for one way, false for the other

// =============================================================
//  VOLTAGE / CURRENT (POWER MONITOR) SETUP
// =============================================================
Adafruit_ADS1115 ads;         // ADS1115 for current (SCT-013)
const int VOLT_PIN = A0;      // ZMPT101B voltage sensor

// Calibration constants — tune to match your multimeter / Kill-A-Watt
float VOLTAGE_CAL  = 1.15;      // Adjust so V reads ~120 V at full power
float CURRENT_CAL  = 0.01028;   // Adjust so I matches a known load
const int V_OFFSET = 511;       // Center point of the AC wave (usually ~512)

// Non-blocking power read interval
unsigned long previousPowerMillis = 0;
const long powerInterval = 1000; // ms between power readings

// RMS sampling window length (200 ms ≈ 12 full cycles of 60 Hz)
const unsigned long rmsSampleWindowMs = 200;

// Latest power readings — updated every powerInterval ms,
// printed every tempInterval ms alongside thermocouple data.
double lastVrms    = 0.0;
double lastIrms    = 0.0;
double lastWattage = 0.0;

// =============================================================
//  HELPERS: keep the pump stepping during blocking operations
// =============================================================

// Single pump-service call — fire one pulse if it's time
inline void servicePump() {
  unsigned long now = micros();
  if (now - previousStepMicros >= stepIntervalMicros) {
    previousStepMicros = now;
    digitalWrite(PUL_PIN, LOW);
    delayMicroseconds(50);
    digitalWrite(PUL_PIN, HIGH);
  }
}

// Keep the pump stepping for a given number of milliseconds
void pumpDuringDelay(unsigned long waitMillis) {
  unsigned long startWait = millis();
  while (millis() - startWait < waitMillis) {
    servicePump();
  }
}

// =============================================================
//  SETUP
// =============================================================
void setup() {
  Serial.begin(115200);
  while (!Serial) delay(1);

  // --- Stepper pins ---
  pinMode(PUL_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  pinMode(ENA_PIN, OUTPUT);

  // Common-anode logic
  digitalWrite(PUL_PIN, HIGH); // Default state HIGH (no pulse)
  digitalWrite(ENA_PIN, HIGH); // HIGH = no current = driver ENABLED

  if (motorDirection) {
    digitalWrite(DIR_PIN, LOW);
  } else {
    digitalWrite(DIR_PIN, HIGH);
  }

  // --- Thermocouples ---
  for (int i = 0; i < numSensors; i++) {
    if (!thermocouples[i]->begin()) {
      Serial.print("Error initializing Thermocouple ");
      Serial.println(i + 1);
      while (1) delay(10);
    }
  }

  // --- ADS1115 (current sensor) ---
  ads.setGain(GAIN_TWO); // ±2.048 V range (good for SCT-013 with burden)
  if (!ads.begin()) {
    Serial.println("Failed to initialize ADS1115. Check wiring!");
    while (1);
  }
  ads.setDataRate(RATE_ADS1115_860SPS);

  // --- CSV header ---
  Serial.println("Time_ms,TC1_C,TC2_C,TC3_C,TC4_C,TC5_C,V_RMS,I_RMS,Wattage");
  pumpDuringDelay(500); // Let MAX31855 stabilize (pump keeps flowing)
}

// =============================================================
//  MAIN LOOP
// =============================================================
void loop() {
  unsigned long currentMillis = millis();
  unsigned long currentMicros = micros();

  // -----------------------------------------------------------
  // 1. NON-BLOCKING STEPPER MOTOR CONTROL
  // -----------------------------------------------------------
  if (currentMicros - previousStepMicros >= stepIntervalMicros) {
    previousStepMicros = currentMicros;
    digitalWrite(PUL_PIN, LOW);
    delayMicroseconds(50);
    digitalWrite(PUL_PIN, HIGH);
  }

  // -----------------------------------------------------------
  // 2. NON-BLOCKING TEMPERATURE READINGS
  // -----------------------------------------------------------
  if (currentMillis - previousTempMillis >= tempInterval) {
    previousTempMillis = currentMillis;

    Serial.print(currentMillis);
    Serial.print(",");

    for (int i = 0; i < numSensors; i++) {
      double temp = thermocouples[i]->readCelsius();

      // Retry logic — pump keeps flowing during retries
      byte retries = 0;
      while (isnan(temp) && retries < 3) {
        pumpDuringDelay(50); // wait 50 ms while still pulsing the pump
        temp = thermocouples[i]->readCelsius();
        retries++;
      }

      if (isnan(temp)) {
        Serial.print("NaN");
      } else {
        Serial.print(temp);
      }

      Serial.print(",");
    }

    // Print latest power values (or 0 if we haven't sampled yet)
    // These will update on their own schedule via the power block below
    Serial.print(lastVrms);
    Serial.print(",");
    Serial.print(lastIrms);
    Serial.print(",");
    Serial.println(lastWattage);
  }

  // -----------------------------------------------------------
  // 3. NON-BLOCKING POWER (V / I) READINGS
  //    The RMS sampling window runs for rmsSampleWindowMs, but
  //    the pump is serviced inside the loop so flow never stops.
  // -----------------------------------------------------------
  if (currentMillis - previousPowerMillis >= powerInterval) {
    previousPowerMillis = currentMillis;

    long v_sq_sum = 0;
    long i_sq_sum = 0;
    long sample_count = 0;

    unsigned long sampleStart = millis();

    while (millis() - sampleStart < rmsSampleWindowMs) {
      // --- keep the pump alive inside the sampling window ---
      servicePump();

      // Voltage sample (ZMPT101B — fast analogRead)
      int16_t v_raw = analogRead(VOLT_PIN) - V_OFFSET;

      // Current sample (SCT-013 via ADS1115 differential)
      int16_t i_raw = ads.readADC_Differential_0_1();

      v_sq_sum += (long)v_raw * v_raw;
      i_sq_sum += (long)i_raw * i_raw;
      sample_count++;
    }

    // Calculate true RMS
    double v_rms = sqrt((double)v_sq_sum / sample_count) * VOLTAGE_CAL;
    double i_rms = sqrt((double)i_sq_sum / sample_count) * CURRENT_CAL;
    double wattage = v_rms * i_rms;

    if (v_rms < 5.0) {
      v_rms = 0.0;
      wattage = 0.0;
    }

    // Cache for printing alongside thermocouple data
    lastVrms    = v_rms;
    lastIrms    = i_rms;
    lastWattage = wattage;
  }
}
