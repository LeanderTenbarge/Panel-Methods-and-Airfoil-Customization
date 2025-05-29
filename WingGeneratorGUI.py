#Multi element wing GUI file: Leander Tenbarge
# Notes:
# This contains a multi element wing creation program that utilizes class shape transformations (Cst), NACA 4 digit camber functions and other things
# to serve as the basis for optimization. While the actual code to the main project will be done in a separate TUI program with different 
# functions (.geo file creation,GMSH scripting, and Openfoam scripting ), this will serve as a visual tool for presenting in FSAE design and to 
# Serve as a tool for future design teams. 

# Loading in the Packages:
import numpy as np
from scipy.special import comb
import tkinter as tk
from matplotlib import pyplot


# Defining the CST foils Class (camber functions, unit scaling, and rotation):
class Foils: 
    # Initialize:
    def __init__(self, p1 , p2 , p3 , p4 ,p5 ,p6 ,p7 ,p8 ,p9 ,p10 ,p11 ,p12,m,p):
        self.parameters_upper = [p1,p2,p3,p4,p5,p6]
        self.parameters_lower = [p7,p8,p9,p10,p11,p12]
        self.x = np.flip(np.geomspace(0.000001,1,40))
        self.y = []
        self.change_thickness()
        self.x_max = np.max(self.x)
        self.y_max = np.max(self.y)
        


    # Apply the CST with the Assumed guesses:
    def change_thickness(self):

        # Base x distribution
        x_dist = np.flip(np.geomspace(0.000001, 1, 40))
        # Parameters
        n = 0.5
        m = 1
        # Upper surface
        y_upper = [((x**n)*(1 - x)**m) * sum(self.parameters_upper[i] * comb(len(self.parameters_upper) - 1, i) * x**i * (1 - x)**(len(self.parameters_upper) - 1 - i) for i in range(len(self.parameters_upper))) for x in x_dist]
        # Lower surface
        y_lower = [-((x**n)*(1 - x)**m) * sum(self.parameters_lower[i] * comb(len(self.parameters_lower) - 1, i) * x**i * (1 - x)**(len(self.parameters_lower) - 1 - i) for i in range(len(self.parameters_lower))) for x in x_dist]
        # Flip lower
        x_lower = np.flip(x_dist[1:])
        y_lower = np.flip(y_lower[1:])
        # Final surfaces
        self.x = np.concatenate((x_dist, x_lower))
        self.y = np.concatenate((y_upper, y_lower))
        self.x = np.append(self.x, self.x[0])
        self.y = np.append(self.y, self.y[0])
        return self.x, self.y
    
    # Utilize the NACA camber function to apply a camber to the airfoil:
    def camber(self,m,p):

        # Enforces the NACA Camber Function:
        for i in range(len(self.x)):
            if self.x[i] > p:
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
        
    
# Defining the Wing instance class, This allows for the calculation of the wing coordinates separately:
class Wing:
    def __init__(self,clr,a,x_off,y_off,x_coords,y_coords):

        if clr == 0:  
            self.valid = False
            self.x = np.array([])
            self.y = np.array([])
            self.area = self.cx = self.cy = self.ixx = self.iyy = 0
            return
        
        self.valid = True 
        self.parameters = [clr,a,x_off,y_off] 
        self.x = x_coords
        self.y = y_coords
        self.rotate()
        self.scale()
        self.x_max = np.max(self.x)
        self.y_max = np.max(self.y)
        self.area = 0
        self.cx = 0
        self.cy = 0 
        self.ixx = 0 
        self.iyy = 0
        self.find_area_parameters()



    #Defining the rotation matrix and the scaling factor:
    def rotate(self):
        r_matrix = np.array([[np.cos(self.parameters[1]) , -np.sin(self.parameters[1])],[np.sin(self.parameters[1]) , np.cos(self.parameters[1])]])
        points = [
        (
            x * r_matrix[0][0] + y * r_matrix[0][1],  
            x * r_matrix[1][0] + y * r_matrix[1][1]  
        )
        for x, y in zip(self.x, self.y)
        ]
        self.x,self.y = zip(*points)
        self.x = np.array(self.x)
        self.y = np.array(self.y)
        
    def scale(self):
        scale = self.parameters[0]
        x_off = self.parameters[2]
        y_off = self.parameters[3]
    
        self.x = self.x * scale + x_off
        self.y = self.y * scale + y_off

    def find_area_parameters(self):
        for i in range(len(self.x)-1):
            self.area += (self.x[i]*self.y[i+1]-self.x[i+1]*self.y[i])
        self.area = 0.5* abs(self.area)

        for i in range(len(self.x)-1):
            self.cx += (self.x[i]+self.x[i+1])*(self.x[i]*self.y[i+1]-self.x[i+1]*self.y[i])
            self.cy += (self.y[i]+self.y[i+1])*(self.x[i]*self.y[i+1]-self.x[i+1]*self.y[i])
        self.cx = (1/(6*self.area))*self.cx
        self.cy = (1/(6*self.area))*self.cy

        for i in range(len(self.x)-1):
            self.ixx += (self.y[i]**2 +self.y[i]*self.y[i+1]+self.y[i+1]**2)*(self.x[i]*self.y[i+1]-self.x[i+1]*self.y[i]) 
            self.iyy += (self.x[i]**2 +self.x[i]*self.x[i+1]+self.x[i+1]**2)*(self.x[i]*self.y[i+1]-self.x[i+1]*self.y[i]) 

        self.ixx = (1/12)*self.ixx - self.area * self.cy**2
        self.iyy = (1/12)*self.iyy - self.area * self.cx**2

# Defining the Run functions for all 3 CST implementations and the Multi element wing File:


#Little function to reset the entry for error validation
def reset_entry(entry):
        entry.delete(0, tk.END)
        entry.config(bg="white")

# Creating the Gui instance and contruction the Frames:
class CST_frame:
    def __init__(self,row,column,title):
        self.row = row
        self.column = column
        self.title = title
        self.frame = tk.Frame(root,bd = 1,relief = 'solid')
        self.frame.grid(row = self.row,column = self.column,sticky='nsew')
        self.column_pos = [0,0,0,0,0,0,1,1,1,1,1,1,1,1]
        self.row_pos = [1,2,3,4,5,6,1,2,3,4,5,6,7,8]
        self.default_values = ['0.17','0.16','0.14','0.17','0.11','0.18','0.17','0.16','0.14','0.17','0.11','0.18','0.00','0.00']
        self.widgets = []
        self.create_entryboxes()


    def create_entryboxes(self):
        # Creates the Labels for the CST Coefficients
        tk.Label(self.frame, text="Upper Coefficients").grid(row=0, column=0, padx=2, pady=2)
        tk.Label(self.frame, text="Lower Coefficients").grid(row=0, column=1, padx=2, pady=2)
        tk.Label(self.frame, text=self.title).grid(row=0, column=2, padx=2, pady=2)

        # Creates the Entry boxes for the CST Coefficients
        for i in range(len(self.row_pos)):
            object = tk.Entry(self.frame)
            object.grid(row=self.row_pos[i],column=self.column_pos[i], padx = 2,pady = 2)
            object.insert(0,self.default_values[i])
            self.widgets.append(object)

        # Creates the labels for the Individual Camber functions
        tk.Label(self.frame, text="M (maximum camber)").grid(row=7, column=0, padx=2, pady=2)
        tk.Label(self.frame, text="P (maximum location)",).grid(row=8, column=0, padx=2, pady=2)

        # Creates the button
        self.button = tk.Button(self.frame,text = "Calculate",command = self.run).grid(row=8,column=2, padx = 2,pady = 2)
        self.canvas_widget = tk.Canvas(self.frame, bg='white')
        self.canvas_widget.grid(row=1, column=2, rowspan=7, padx=2, pady=2, sticky='nsew')

    # Run command
    def run(self):

        # Creates the Input parameter array
        self.parameters = []

        # Gets the Values from the text entries
        for i in range(len(self.widgets)):
            try:
                self.parameters.append(float(self.widgets[i].get()))

            # Insures that errors are handled properly
            except ValueError:
                self.widgets[i].delete("0", tk.END)
                self.widgets[i].insert(0,"Error")
                self.widgets[i].config(bg="lightcoral")
                self.widgets[i].after(500, reset_entry, self.widgets[i])
              
        # Actually Executing the code for the CST foil
        foil = Foils(*self.parameters)
        self.x,self.y = foil.camber(self.parameters[-2],self.parameters[-1])

        #plots the Canvas
        plotx = self.x * 250 + 15
        ploty = -1*self.y * 250 + 115
        self.canvas_widget.delete("all")  
        for i in range(len(self.y)):
            self.canvas_widget.create_oval(plotx[i]-2,ploty[i]-2,plotx[i]+2,ploty[i]+2,fill='blue',outline ='blue') 


class canvas_frame:
    def __init__(self,row,column):
        self.row = row
        self.column = column
        self.frame = tk.Frame(root,bd = 1,relief = 'solid')
        self.frame.grid(row = self.row,column = self.column,sticky='nsew')
        self.canvas = tk.Canvas(self.frame, height=275, width=525, bg='white')
        self.canvas.grid(row=0, column=0, columnspan=1, padx=5, pady=5)

class multi_element_parameters:
    def __init__(self,row,column,title):
        self.row = row
        self.column = column
        self.title = title
        self.frame = tk.Frame(root,bd = 1,relief = 'solid')
        self.frame.grid(row = self.row,column = self.column,sticky='nsew')
        self.column_pos = [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2,2,2,2,2,2]
        self.row_pos    = [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2,3,4,5,6,7,8] 
        self.default_values = ['0.60','0.25','0.15','0.05','0.05','0.00','0.10','0.20','0.00','0.00','0.00','0.00','0.00','0.00','0.00','0.00','0.00','0.00','0.00']
        self.names = ['CLR1:','CLR2:','CLR3:','Unit Slotgap:','Unit Overlap:','AOA 1:','AOA 2:','AOA 3:', 'Height/Width','Area 1:', 'Area 2:','Area 3:','Ixx,c 1','Iyy,c 1','Ixx,c 2','Iyy,c 2','Ixx,c 3','Iyy,c 3']
        self.widgets = []
        for i in range(len(self.names)):
            tk.Label(self.frame, text=self.names[i]).grid(row=self.row_pos[i], column=self.column_pos[i], padx=2, pady=2)
            object = tk.Entry(self.frame)
            object.grid(row=self.row_pos[i],column=self.column_pos[i]+1, padx = 2,pady = 2)
            object.insert(0,self.default_values[i])
            self.widgets.append(object)

class textbox:
    def __init__(self,row,column):
        self.row = row
        self.column = column
        self.frame = tk.Frame(root)
        self.frame = tk.Frame(root,bd = 1,relief = 'solid')
        self.button = tk.Button(self.frame,text = "Calculate",command=self.run).grid(row=0,column=0, padx = 2,pady = 2)
        self.frame.grid(row=self.row, column=self.column, columnspan = 1, sticky="n")
        tk.Label(self.frame, text="Primary foil").grid(row=1, column=0)
        self.primary_output = tk.Text(self.frame, width=50, height=2)
        self.primary_output.grid(row=2, column=0, columnspan=1, padx=5, pady=0)
        tk.Label(self.frame, text="Secondary foil").grid(row=3, column=0)
        self.secondary_output = tk.Text(self.frame, width=50, height=2)
        self.secondary_output.grid(row=4, column=0, columnspan=1, padx=5, pady=0)
        tk.Label(self.frame, text="Tertiary foil").grid(row=5, column=0)
        self.tertiary_output = tk.Text(self.frame, width=50, height=2)
        self.tertiary_output.grid(row=6, column=0, columnspan=1, padx=5, pady=0)
   
    def run(self):

        # Creates the Input parameter array
        self.parameters = []

        # Gets the Values from the text entries
        for i in range(len(panel4.widgets)):
            try:
                self.parameters.append(float(panel4.widgets[i].get()))
            
            # Insures that errors are handled properly
            except ValueError:
                panel4.widgets[i].delete("0", tk.END)
                panel4.widgets[i].insert(0,"Error")
                panel4.widgets[i].config(bg="lightcoral")
                panel4.widgets[i].after(500, reset_entry, panel4.widgets[i])

        w1 = Wing(float(panel4.widgets[0].get()),float(panel4.widgets[5].get())*np.pi/180,0,0,panel1.x,panel1.y)
        w2 = Wing(float(panel4.widgets[1].get()),float(panel4.widgets[6].get())*np.pi/180,w1.x[0]-float(panel4.widgets[4].get()),w1.y[0]+float(panel4.widgets[3].get()),panel2.x,panel2.y)
        if panel4.widgets[2].get() != 0:
            w3 = Wing(float(panel4.widgets[2].get()),float(panel4.widgets[7].get())*np.pi/180,w2.x[0]-float(panel4.widgets[4].get()),w2.y[0]+float(panel4.widgets[3].get()),panel3.x,panel3.y)
        else:
            w3 = Wing(0, 0, 0, 0, panel3.x, panel3.y)
            w3.x = np.array([0])
            w3.y = np.array([0])
            
    
        panel4.widgets[9].delete("0", tk.END)
        panel4.widgets[10].delete("0", tk.END)
        panel4.widgets[11].delete("0", tk.END)
        panel4.widgets[12].delete("0", tk.END)
        panel4.widgets[13].delete("0", tk.END)
        panel4.widgets[14].delete("0", tk.END)
        panel4.widgets[15].delete("0", tk.END)
        panel4.widgets[16].delete("0", tk.END)
        panel4.widgets[17].delete("0", tk.END)
        panel4.widgets[9].insert('end', str(w1.area))
        panel4.widgets[10].insert('end',str(w2.area))
        panel4.widgets[11].insert('end',str(w3.area))
        panel4.widgets[12].insert('end',str(w1.ixx))
        panel4.widgets[13].insert('end',str(w1.iyy))
        panel4.widgets[14].insert('end',str(w2.ixx))
        panel4.widgets[15].insert('end',str(w2.iyy))
        panel4.widgets[16].insert('end',str(w3.ixx))
        panel4.widgets[17].insert('end',str(w3.iyy))
        

        panel5.primary_output.delete('1.0',tk.END)
        panel5.secondary_output.delete('1.0',tk.END)
        panel5.tertiary_output.delete('1.0',tk.END)
        panel_canvas.canvas.delete("all")  
        if w1.valid:
            panel4.widgets[9].delete("0", tk.END)
            panel4.widgets[9].insert('end', str(w1.y[0]/w1.x[0]))
            for x, y in zip(w1.x, w1.y):
                panel5.primary_output.insert(tk.END, f"{x:.7f}\t{y:.7f}\n")
                canvas_x = 525 * x+25
                canvas_y = -275 * y + 135  # Apply offset after flipping
                panel_canvas.canvas.create_oval(canvas_x-2, canvas_y-2, canvas_x+2, canvas_y+2, fill='blue', outline='blue')

        if w2.valid:
            panel4.widgets[9].delete("0", tk.END)
            panel4.widgets[9].insert('end', str(w2.y[0]/w2.x[0]))
            for x, y in zip(w2.x, w2.y):
                panel5.secondary_output.insert(tk.END, f"{x:.7f}\t{y:.7f}\n")
                canvas_x = 525 * x
                canvas_y = -275 * y + 135
                panel_canvas.canvas.create_oval(canvas_x-2, canvas_y-2, canvas_x+2, canvas_y+2, fill='blue', outline='blue')

        if w3.valid:
            panel4.widgets[8].delete("0", tk.END)
            panel4.widgets[8].insert('end', str(w3.y[0]/w3.x[0]))
            for x, y in zip(w3.x, w3.y):
                panel5.tertiary_output.insert(tk.END, f"{x:.7f}\t{y:.7f}\n") 
                canvas_x = 525 * x
                canvas_y = -275 * y + 135
                panel_canvas.canvas.create_oval(canvas_x-2, canvas_y-2, canvas_x+2, canvas_y+2, fill='blue', outline='blue')
        
        


#Creates The Gui
root=tk.Tk()
root.title("Multi Element Wing Generator")
root.grid_columnconfigure(0, weight=1)
root.grid_columnconfigure(1, weight=1)

panel1 = CST_frame(0,0,"Airfoil 1")
panel2 = CST_frame(1,0,"Airfoil 2")
panel3 = CST_frame(2,0,"Airfoil 3")
panel4 = multi_element_parameters(0,1,"Parameters")
panel5 = textbox(2,1)
panel_canvas = canvas_frame(1,1)
root.mainloop()











