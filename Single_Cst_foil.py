import numpy as np
from matplotlib import pyplot
from scipy.special import comb
import tkinter as tk

#Class Shape Transformation for a single airfoil with unique upper and lower coordinates.
class Cst_foils:  
    def __init__(self, p1 , p2 , p3 , p4 ,p5 ,p6 ,p7 ,p8 ,p9 ,p10 ,p11 ,p12):
        self.parameters_upper = [p1,p2,p3,p4,p5,p6]
        self.parameters_lower = [p7,p8,p9,p10,p11,p12]
        self.x = np.flip(np.geomspace(0.000001,1,150))
        self.y = []

    def change_thickness(self):
        n = 0.5
        m = 1
        # Upper surface
        y_upper = [((x**n)*(1 - x)**m) * sum(self.parameters_upper[i] * comb(len(self.parameters_upper) - 1, i) * x**i * (1 - x)**(len(self.parameters_upper) - 1 - i) for i in range(len(self.parameters_upper))) for x in self.x]
        # Lower surface
        y_lower = [-((x**n)*(1 - x)**m) * sum(self.parameters_lower[i] * comb(len(self.parameters_lower) - 1, i) * x**i * (1 - x)**(len(self.parameters_lower) - 1 - i) for i in range(len(self.parameters_lower))) for x in self.x]
        # Flip lower surface for trailing edge continuity
        x_lower = np.flip(self.x[1:])          
        y_lower = np.flip(y_lower[1:])
        # Combine upper and lower
        self.x = np.concatenate((self.x, x_lower))
        self.y = np.concatenate((y_upper, y_lower))
        return self.x, self.y
    
    def camber(self,m,p):
        for i in range(len(self.x)):
            if self.x[i] >= p:
                try:
                    self.y[i] += (m/p**2)*(2*p*self.x[i]-self.x[i]**2)
                except:
                    pass
            else:
                try:
                    self.y[i] += (m/(1-p)**2)*((1-2*p)+2*p*self.x[i]-(self.x[i]**2))
                except:
                    pass
        return self.x,self.y
    
def run():
    p1 = float(ti1.get())
    p2 = float(ti2.get())
    p3 = float(ti3.get())
    p4 = float(ti4.get())
    p5 = float(ti5.get())
    p6 = float(ti6.get())
    p7 = float(ti7.get())
    p8 = float(ti8.get())
    p9 = float(ti9.get())
    p10 = float(ti10.get())
    p11 = float(ti11.get())
    p12 = float(ti12.get())
    p13 = float(ti13.get())
    p14 = float(ti14.get())
    Alpha = Cst_foils(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
    ax,ay = Alpha.change_thickness()
    ax,ay = Alpha.camber(p13,p14)
    plotx = ax * 640 + 10
    ploty = -1*ay * 500 + 150
    canvas.delete("all")  
    for i in range(len(ay)):
        canvas.create_oval(plotx[i]-2,ploty[i]-2,plotx[i]+2,ploty[i]+2,fill='blue',outline ='blue')

    # Clear previous output
    dat_output.delete("1.0", tk.END)

    # Write new coordinate data to the Coordinate Display
    for x, y in zip(ax, ay):
        dat_output.insert(tk.END, f"{x:.7f}\t{y:.7f}\n")

    #Fixes the trailing edge issue 
    dat_output.insert(tk.END, f"{ax[0]:.7f}\t{ay[0]:.7f}\n")


#Creates The Gui
root=tk.Tk()
root.title("Single Airfoil CST implementation")
# Formats the Text inputs in one frame

text_frame = tk.Frame(root)
text_frame.grid(row=0, column=0, padx=10, pady=10, sticky='n')

#Labels and inputs to the first panel
tk.Label(text_frame, text="Upper Coefficients").grid(row=0, column=0, padx=5, pady=5)
tk.Label(text_frame, text="Lower Coefficients").grid(row=0, column=1, padx=5, pady=5)
ti1 = tk.Entry(text_frame)
ti1.grid(row=1,column=0, padx = 5,pady = 5)
ti1.insert(0,"0.17")
ti2 = tk.Entry(text_frame)
ti2.grid(row=2,column=0, padx = 5,pady = 5)
ti2.insert(0,"0.16")
ti3 = tk.Entry(text_frame)
ti3.grid(row=3,column=0, padx = 5,pady = 5)
ti3.insert(0,"0.14")
ti4 = tk.Entry(text_frame)
ti4.grid(row=4,column=0, padx = 5,pady = 5)
ti4.insert(0,"0.17")
ti5 = tk.Entry(text_frame)
ti5.grid(row=5,column=0, padx = 5,pady = 5)
ti5.insert(0,"0.11")
ti6 = tk.Entry(text_frame)
ti6.grid(row=6,column=0, padx = 5,pady = 5)
ti6.insert(0,"0.18")
ti7 = tk.Entry(text_frame)
ti7.grid(row=1,column=1, padx = 5,pady = 5)
ti7.insert(0,"0.17")
ti8 = tk.Entry(text_frame)
ti8.grid(row=2,column=1, padx = 5,pady = 5)
ti8.insert(0,"0.16")
ti9 = tk.Entry(text_frame)
ti9.grid(row=3,column=1, padx = 5,pady = 5)
ti9.insert(0,"0.14")
ti10 = tk.Entry(text_frame)
ti10.grid(row=4,column=1, padx = 5,pady = 5)
ti10.insert(0,"0.17")
ti11 = tk.Entry(text_frame)
ti11.grid(row=5,column=1, padx = 5,pady = 5)
ti11.insert(0,".11")
ti12 = tk.Entry(text_frame)
ti12.grid(row=6,column=1, padx = 5,pady = 5)
ti12.insert(0,".18")
tk.Label(text_frame, text="Camber Functions").grid(row=7, column=0, padx=5, pady=5)
tk.Label(text_frame, text="M (maximum camber)",).grid(row=8, column=0, padx=5, pady=5)
ti13 = tk.Entry(text_frame)
ti13.grid(row=8,column=1, padx = 5,pady = 5)
ti13.insert(0,"0")
tk.Label(text_frame, text="P  (maximum location)").grid(row=9, column=0, padx=5, pady=5)
ti14 = tk.Entry(text_frame)
ti14.grid(row=9,column=1, padx = 5,pady = 5)
ti14.insert(0,"0")
b1 = tk.Button(text_frame,text = "Calculate",command = run,).grid(row=10,column=0, padx = 5,pady = 5)

# Creating the Frame for the Canvas
canvas_frame = tk.Frame(root)
canvas_frame.grid(row=0, column=4, columnspan=3, sticky="n")
tk.Label(canvas_frame, text="Airfoil Plot").grid(row=0, column=0, columnspan=2, pady=(0, 2))
canvas = tk.Canvas(canvas_frame, height=390, width=650, bg='white')
canvas.grid(row=1, column=0, columnspan=2, padx=5, pady=0)

# Creating the Panel For the Dat File Scrolling
dat_frame = tk.Frame(root)
dat_frame.grid(row=0, column=7, columnspan=3, sticky="n")
tk.Label(dat_frame, text="Airfoil Coordinates (text delimited)").grid(row=0, column=1, sticky="w")
dat_output = tk.Text(dat_frame, width=50, height=30)
dat_output.grid(row=1, column=0, columnspan=2, padx=5, pady=0)
root.mainloop()



