import serial
import csv
import time

# --- SERIAL CONFIGURATION ---
# 1. Find your Arduino's COM port (check Device Manager)
# 2. Set the baud rate to match the Arduino code (115200)
SERIAL_PORT = 'COM4'  # Change this to your Arduino's port
BAUD_RATE = 115200

# --- FILE CONFIGURATION ---
# The name of the CSV file to create
CSV_FILENAME = 'reactor_data.csv'

# Column headers (must match the Arduino output format)
HEADERS = ['Time_ms', 'TC1_C', 'TC2_C', 'TC3_C', 'TC4_C', 'TC5_C']

def main():
    print(f"Attempting to connect to {SERIAL_PORT} at {BAUD_RATE} baud...")
    
    try:
        # Initialize serial connection
        ser = serial.Serial(SERIAL_PORT, BAUD_RATE, timeout=1)
        time.sleep(2) # Wait for the connection to establish
        
        # Open the CSV file for writing
        with open(CSV_FILENAME, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write the header row
            writer.writerow(HEADERS)
            print(f"Successfully opened {CSV_FILENAME}. Starting data collection...")
            print("Press Ctrl+C to stop logging.")

            while True:
                # Read a line from the serial port
                line = ser.readline()
                
                # Decode bytes to string and strip whitespace
                try:
                    line_str = line.decode('utf-8').strip()
                except UnicodeDecodeError:
                    continue # Skip if data is corrupted

                # Check if the line is empty
                if line_str:
                    # Split the string by comma to get individual values
                    values = line_str.split(',')
                    
                    # Write the values to the CSV file
                    writer.writerow(values)
                    
                    # Also print to console for real-time feedback
                    print(line_str)

    except serial.SerialException as e:
        print(f"Error: {e}")
        print("Please check if the Arduino is connected and the COM port is correct.")
    except KeyboardInterrupt:
        print("\nData collection stopped by user.")
    finally:
        if 'ser' in locals() and ser.is_open:
            ser.close()
            print("Serial connection closed.")

if __name__ == "__main__":
    main()