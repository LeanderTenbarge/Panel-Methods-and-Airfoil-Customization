#Leander Tenbarge
#GUI for a CST parameter curve fit

# Imports
import tkinter as tk
import numpy as np
from scipy.special import comb
from scipy.optimize import curve_fit

#The Data Parsing Function
def parse_airfoil_data(text_widget):
    raw = text_widget.get("1.0", tk.END).strip()
    coords = []

    for line in raw.splitlines():
        try:
            x_str, y_str = line.strip().split()
            x, y = float(x_str), float(y_str)
            coords.append((x, y))
        except ValueError:
            continue  # skip bad lines

    coords = np.array(coords)
    #split at the x=0 point or midpoint of array
    mid = len(coords) // 2
    upper = coords[:mid]
    lower = coords[mid:]

    # Sort x increasing for fitting
    upper = upper[np.argsort(upper[:, 0])]
    lower = lower[np.argsort(lower[:, 0])]

    return upper[:, 0], upper[:, 1], lower[:, 0], lower[:, 1]

# Class Function
def class_func(x, N=0.5, M=1.0):
    return (x**N) * ((1 - x)**M)

# Bernstein Basis function
def bernstein(i, n, x):
    return comb(n, i) * x**i * (1 - x)**(n - i)

# CST model
def cst_model(x, *A):
    n = len(A) - 1
    shape = sum(A[i] * bernstein(i, n, x) for i in range(n + 1))
    return class_func(x) * shape


# Setting up the actual function
def run():
    U1.delete(0,tk.END)
    U2.delete(0,tk.END)
    U3.delete(0,tk.END)
    U4.delete(0,tk.END)
    U5.delete(0,tk.END)
    U6.delete(0,tk.END)
    L1.delete(0,tk.END)
    L2.delete(0,tk.END)
    L3.delete(0,tk.END)
    L4.delete(0,tk.END)
    L5.delete(0,tk.END)
    L6.delete(0,tk.END)
    x_u, y_u, x_l, y_l = parse_airfoil_data(text_entry)
    n = 5
    initial_guess = [0.1] * (n + 1)
    # Fit upper and lower surfaces separately
    c_upper, _ = curve_fit(cst_model, x_u, y_u, p0=initial_guess)
    c_lower, _ = curve_fit(cst_model, x_l, y_l, p0=initial_guess) 
    c_lower = c_lower[:]*-1
    U1.insert(0,f"{(c_upper[0]):.2f}")
    U2.insert(0,f"{(c_upper[1]):.2f}")
    U3.insert(0,f"{(c_upper[2]):.2f}")
    U4.insert(0,f"{(c_upper[3]):.2f}")
    U5.insert(0,f"{(c_upper[4]):.2f}")
    U6.insert(0,f"{(c_upper[5]):.2f}")
    L1.insert(0,f"{(c_lower[0]):.2f}")
    L2.insert(0,f"{(c_lower[1]):.2f}")
    L3.insert(0,f"{(c_lower[2]):.2f}")
    L4.insert(0,f"{(c_lower[3]):.2f}")
    L5.insert(0,f"{(c_lower[4]):.2f}")
    L6.insert(0,f"{(c_lower[5]):.2f}")
  

# Setting up the Gui
root = tk.Tk()
root.title("CST parameter solver")

# Setting up the frame for the Text input 
text_input = tk.Frame(root)
text_input.grid(row=0, column=0, columnspan = 3,padx=10, pady=10, sticky='nsew')
text_label = tk.Label(text_input,text = "Selig DAT format input").grid(row = 0,column = 0 ,padx = 5, pady = 5,sticky = 'nsew' )
text_entry = tk.Text(text_input, width=50, height=30)
text_entry.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky='nsew')


# Setting up the Calculation Button
button = tk.Button(text_input,text = "Calculate",command = run).grid(row=0, column=1)

# Setting up the Coefficients Frame
c_frame = tk.Frame(root)
c_frame.grid(row = 0,column = 3,columnspan = 3,padx = 10, pady = 10 , sticky = 'nsew')
c_label_u = tk.Label(c_frame,text = 'Upper coefficients').grid(row = 0 ,column = 0 , padx = 5, pady =5)
c_label_l = tk.Label(c_frame,text = 'Lower coefficients').grid(row = 0 ,column = 1 , padx = 5, pady =5)

#Creating the input labels 
U1 = tk.Entry(c_frame)
U2 = tk.Entry(c_frame)
U3 = tk.Entry(c_frame)
U4 = tk.Entry(c_frame)
U5 = tk.Entry(c_frame)
U6 = tk.Entry(c_frame)
L1 = tk.Entry(c_frame)
L2 = tk.Entry(c_frame)
L3 = tk.Entry(c_frame)
L4 = tk.Entry(c_frame)
L5 = tk.Entry(c_frame)
L6 = tk.Entry(c_frame)
U1.grid(row = 1,column = 0,padx = 5, pady = 5)
U2.grid(row = 2,column = 0,padx = 5, pady = 5)
U3.grid(row = 3,column = 0,padx = 5, pady = 5)
U4.grid(row = 4,column = 0,padx = 5, pady = 5)
U5.grid(row = 5,column = 0,padx = 5, pady = 5)
U6.grid(row = 6,column = 0,padx = 5, pady = 5)
L1.grid(row = 1,column = 1,padx = 5, pady = 5)
L2.grid(row = 2,column = 1,padx = 5, pady = 5)
L3.grid(row = 3,column = 1,padx = 5, pady = 5)
L4.grid(row = 4,column = 1,padx = 5, pady = 5)
L5.grid(row = 5,column = 1,padx = 5, pady = 5)
L6.grid(row = 6,column = 1,padx = 5, pady = 5)
N1 = tk.Text(c_frame,width=20, height=15)
N1.grid(row = 6,column = 0,columnspan = 3, padx = 10 , pady = 10, sticky = 'nsew' )
N1.insert("1.0",'Note: The input must be in the Selig-format .dat file, which can be obtained by visiting http://airfoiltools.com/airfoil, selecting an airfoil, and choosing the "Selig format Dat File" option. Be sure to omit the title header line from the file. Also, please note that there is no input error validation — all values must be correctly formatted and valid numbers.')
root.mainloop()