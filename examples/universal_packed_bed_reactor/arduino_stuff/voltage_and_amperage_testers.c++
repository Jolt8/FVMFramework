#include <Wire.h>
#include <Adafruit_ADS1X15.h>

// --- HARDWARE CONFIG ---
Adafruit_ADS1115 ads;     // ADS1115 for Current (SCT-013)
const int VOLT_PIN = A0;  // ZMPT101B Voltage Sensor

// --- CALIBRATION CONSTANTS ---
// You will tune these to match your multimeter/Kill-A-Watt
float VOLTAGE_CAL = 0.585;  // Adjust so V reads ~120V at full power
float CURRENT_CAL = 0.042;  // Adjust so I matches a known load
const int V_OFFSET = 512;   // The center point of the AC wave (usually 512)

void setup() {
  Serial.begin(115200);
  
  // Initialize ADS1115
  ads.setGain(GAIN_TWO);    // +/- 2.048V range (Good for SCT-013 with burden)
  if (!ads.begin()) {
    Serial.println("Failed to initialize ADS1115. Check wiring!");
    while (1);
  }

  ads.setDataRate(RATE_ADS1115_860SPS)

  Serial.println("Power Monitor Test Started...");
  Serial.println("V_RMS, I_RMS, Wattage");
}

void loop() {
  uint32_t start_time = millis();
  
  long v_sq_sum = 0;
  long i_sq_sum = 0;
  long sample_count = 0;

  // --- SAMPLING WINDOW ---
  // We sample for 50ms (3 full cycles of 60Hz) to get an accurate RMS
  while (millis() - start_time < 50) {
    // Read Voltage (ZMPT101B)
    int16_t v_raw = analogRead(VOLT_PIN) - V_OFFSET;

    Serial.print(analogRead(VOLT_PIN))
    
    // Read Current (SCT-013 Differential)
    // A0 and A1 pins on the ADS1115
    int16_t i_raw = ads.readADC_Differential_0_1();

    v_sq_sum += (long)v_raw * v_raw;
    i_sq_sum += (long)i_raw * i_raw;
    sample_count++;
  }

  // --- CALCULATE TRUE RMS ---
  double v_rms = sqrt(v_sq_sum / sample_count) * VOLTAGE_CAL;
  double i_rms = sqrt(i_sq_sum / sample_count) * CURRENT_CAL;
  
  // Power (P = V * I)
  double wattage = v_rms * i_rms;

  if (v_rms < 5.0) { 
    v_rms = 0.0; 
    wattage = 0.0;
  }

  // Print results
  Serial.print(v_rms);
  Serial.print("V, ");
  Serial.print(i_rms);
  Serial.print("A, ");
  Serial.print(wattage);
  Serial.println("W");

  delay(500); // Wait half a second before next reading
}