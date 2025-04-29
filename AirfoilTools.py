#Leander Tenbarge
# This is a collection of Tools used to manipulate Airfoils
import sys
import numpy as np



# Returns a NumPy array with the respective x and y coordinates (x, y)
# Files are in the SELIG .dat format as found in the UIUC Airfoil database 
def get_coordinates(filename):
    xcoords = []
    ycoords = []

    if isinstance(filename, str):
        with open(filename, 'r') as file:
            next(file)  # skip the first line (airfoil name)
            for line in file:
                if line.strip():  # skip blank lines
                    parts = line.strip().split()  # split on any whitespace
                    if len(parts) >= 2:
                        try:
                            x, y = map(float, parts[:2])
                            xcoords.append(x)
                            ycoords.append(y)
                        except ValueError:
                            print(f"Skipping malformed line: {line}")
                            continue

        return np.array(xcoords), np.array(ycoords)

    else:
        print("The filename provided is not a string.")
        sys.exit(1)
        
