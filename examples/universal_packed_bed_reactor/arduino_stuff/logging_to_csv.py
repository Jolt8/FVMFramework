import serial
import csv
import time
import os
from datetime import datetime

# --- SERIAL CONFIGURATION ---
SERIAL_PORT = 'COM4'  # Change this to your Arduino's port
BAUD_RATE = 115200

# Column headers
HEADERS = ['Time_ms', 'TC1_C', 'TC2_C', 'TC3_C', 'TC4_C', 'TC5_C']

def main():
    # --- GENERATE DYNAMIC FILENAME ---
    # Format: reactor_data_YYYY_MM_DD__HH_MM_SS.csv
    # Added seconds just in case you restart the script twice in one minute!
    now = datetime.now()
    csv_filename = now.strftime("reactor_data_%Y_%m_%d__%H_%M_%S.csv")
    
    print(f"Attempting to connect to {SERIAL_PORT} at {BAUD_RATE} baud...")
    
    try:
        # Initialize serial connection
        ser = serial.Serial(SERIAL_PORT, BAUD_RATE, timeout=1)
        time.sleep(2) # Wait for the Arduino to reset and connect
        
        # Open the CSV file for writing
        with open(csv_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write the header row
            writer.writerow(HEADERS)
            print(f"Successfully opened {csv_filename}.")
            print("--- DATA LOGGING LIVE ---")
            print("Press Ctrl+C to stop logging.")

            # Flush the header to disk immediately
            csvfile.flush()
            os.fsync(csvfile.fileno())

            while True:
                # Read a line from the serial port
                line = ser.readline()
                
                # Decode bytes to string and strip whitespace
                try:
                    line_str = line.decode('utf-8').strip()
                except UnicodeDecodeError:
                    continue # Skip if electrical noise corrupts a byte

                # Ignore empty lines or the Arduino's text header
                if not line_str or "Time" in line_str or "Error" in line_str:
                    continue 

                # Check if the line actually has data
                if line_str:
                    values = line_str.split(',')
                    
                    # Sanity check: Did we get exactly 6 pieces of data?
                    # This prevents half-written lines from ruining your CSV columns
                    if len(values) == len(HEADERS):
                        # 1. Write to the CSV buffer
                        writer.writerow(values)
                        
                        # 2. FORCE WRITE TO HARD DRIVE (The 100% Reliability Hack)
                        csvfile.flush()              # Push from Python to OS
                        os.fsync(csvfile.fileno())   # Push from OS to Physical Disk
                        
                        # 3. Print to console
                        print(line_str)
                    else:
                        print(f"[IGNORED MALFORMED LINE]: {line_str}")

    except serial.SerialException as e:
        print(f"\n[HARDWARE ERROR]: {e}")
        print("Check if the Arduino is connected, or if another program (like Arduino IDE) is using the COM port.")
    except KeyboardInterrupt:
        print(f"\n[SUCCESS]: Data collection stopped cleanly by user. Saved to {csv_filename}")
    finally:
        if 'ser' in locals() and ser.is_open:
            ser.close()
            print("Serial connection safely closed.")

if __name__ == "__main__":
    main()