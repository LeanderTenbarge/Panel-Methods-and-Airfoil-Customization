
#Source Panel formulation
#Leander Tenbarge
#Last updated - 3/16/2025
#Note: This code is written with goal of obtaining pressure coefficent graphs with the intention of comparing with -
# - other panel methods. While the code has the ability to analyse at different Angles of attack, it is imperative to understand - 
# - mathematical formulations behind the form of analysis. The Flow is inviscid, irrotational and nonlifting, we can only - 
# - indirectly gather lift values due to the pressure distribution.

#Importing modules
import sys
import numpy as np
from matplotlib import pyplot , path
import math

# Check for Counter clockwise orientation utilizing the Gauss's Area Formulation: 
def check_direction(x_bound,y_bound):
    if len(x_bound) == len(y_bound):
        x_bound = np.array(x_bound)
        y_bound = np.array(y_bound)
        area = (1/2)*(np.sum(x_bound[0:-1]*y_bound[1:])-sum(x_bound[1:]*y_bound[0:-1]))

        if area > 0:
            x_bound = x_bound[::-1]
            y_bound = y_bound[::-1]
    else:
        print("The X and Y lists are unequal in length")
        sys.exit()


    return x_bound, y_bound

#Process the input geometry, calculate Central points, Calculate angles, Adjust angle outputs to be (+):  
def process_geometry(x_bound,y_bound):
    
    #Initalizing the NP arrays
    x_cont = np.zeros(len(x_bound) - 1 )
    y_cont = np.zeros(len(x_bound) - 1 )
    s_length = np.zeros(len(x_bound) - 1 )
    phi = np.zeros(len(x_bound) - 1 )


    #Calculating the Control points (midpoint) of each panel
    for i in range(len(x_bound) - 1 ):
        x_cont[i] = (x_bound[i]+ x_bound[i+1])/(2)
    for i in range(len(y_bound) - 1):
        y_cont[i] = (y_bound[i] + y_bound[i+1])/(2)

    #Calculating the length of each panel
    for i in range(len(s_length)):
        s_length[i] = np.sqrt((x_bound[i+1]- x_bound[i])**2 + (y_bound[i+1]- y_bound[i])**2)

    #Calculating the Panel angles (RAD) and converting all of the panels to positive
    for i in range(len(phi)):
        dx = x_bound[i+1]-x_bound[i]
        dy = y_bound[i+1]-y_bound[i]
        phi[i] = np.atan2(dy,dx)
        if phi[i] < 0:
            phi[i] = phi[i] + 2*np.pi

    return x_cont , y_cont , s_length , phi

#Given I and J panel find The Geometric integral I_ij & J_ij: Status: 
def sp_geometric_integrals(x_cont,y_cont,x_bound,y_bound,phi,s_length):

    #Initalizing the Geometric integral arrays
    I = np.zeros([len(x_cont),len(x_cont)])
    J = np.zeros([len(x_cont),len(x_cont)])
    
    #looping over the matricies
    for i in range(len(x_cont)):
        for j in range(len(x_cont)):
            if i != j:

                #Solving for the coefficients 
                a = - (x_cont[i] - x_bound[j])*np.cos(phi[j]) - (y_cont[i]-y_bound[j])*np.sin(phi[j])
                b = (x_cont[i]-x_bound[j])**2+(y_cont[i]-y_bound[j])**2
                c_n = np.sin(phi[i]-phi[j])
                d_n = -(x_cont[i] - x_bound[j])* np.sin(phi[i]) + (y_cont[i]-y_bound[j])*np.cos(phi[i])
                c_t = - np.cos(phi[i]-phi[j])
                d_t = (x_cont[i]-x_bound[j])*np.cos(phi[i]) + (y_cont[i] - y_bound[j])*np.sin(phi[j])
                e = np.sqrt(b-a**2)

                #Finding the arrays of geometric integrals
                I[i,j] = (c_n/2)*(np.log((s_length[j]**2 + 2*a*s_length[j] + b)/(b))) + ((d_n-a*c_n)/(e))*(math.atan2((s_length[j]+a),(e)) - math.atan2((a),(e)))
                J[i,j] = (c_t/2)*(np.log((s_length[j]**2 + 2*a*s_length[j] + b)/(b))) + ((d_t-a*c_t)/(e))*(math.atan2((s_length[j]+a),(e)) - math.atan2((a),(e)))

                # Check if the value is NANinf or divide by zero to mitigate chance of numerical error
                if (np.iscomplex(I[i,j]) or np.isnan(I[i,j]) or np.isinf(I[i,j])):      
                    I[i,j] = 0                                                          
                if (np.iscomplex(J[i,j]) or np.isnan(J[i,j]) or np.isinf(J[i,j])):      
                    J[i,j] = 0             

    
    return I , J
       
#Determine the Panel Strengths: Status: 
def find_source_panel_strengths(I , phi , alpha , v_inf):
    
    #Creating a placeholder with the correct values, Contructing the A matrix in x = a^-1*b
    I_calc = I + (np.pi)*(np.eye(len(I)))

    #Constructing the b(i) array 
    beta = np.zeros(len(I))
    for i in range(len(phi)):
        beta[i] = phi[i] + (np.pi/2) - alpha

    #Contructing the B matrix in x = a^-1*b
    b_calc = np.zeros(len(I))
    for i in range(len(I)):
        b_calc[i] = -v_inf * 2 * np.pi * np.cos(beta[i])

    #solving the equation
    panel_strengths = np.zeros(len(I_calc))
    panel_strengths = np.linalg.solve(I_calc,b_calc)

    return beta  , panel_strengths

#Calculate and initialize the Vector field: 
def sp_calc_point_geometric(N_vectors,x_dim,y_dim,x_bound,y_bound,phi ,s_length):

    #Initializing the arrays
    N_panels = int(len(x_bound) - 1 )
    x = np.linspace( x_dim[0] , x_dim[1] , int(((x_dim[1]-x_dim[0])/(y_dim[1]-y_dim[0])) * N_vectors))
    y = np.linspace( y_dim[0] , y_dim[1] , N_vectors)
    z = np.linspace(0 , N_panels-1 , N_panels)
    x_grid , y_grid = np.meshgrid(x,y,indexing='xy')
    Mx = np.zeros((len(y),len(x),len(z)))
    My = np.zeros((len(y),len(x),len(z)))
   

    # Creating the geomtric integral 
    for i in range(len(y)): # X points
        for j in range(len(x)): # Y points
            for k in range(len(z)):  # Iterating over every single panel for each point:
            
                #Calculating the Coeffcients 
                a = - (x_grid[i,j] - x_bound[k])*np.cos(phi[k]) - (y_grid[i,j] - y_bound[k])*np.sin(phi[k])
                b = (x_grid[i,j] - x_bound[k])**2 + (y_grid[i,j] - y_bound[k])**2
                c_x = -np.cos(phi[k])
                c_y = -np.sin(phi[k])
                d_x = (x_grid[i,j] - x_bound[k])
                d_y = ((y_grid[i,j] - y_bound[k]))
                e = np.sqrt(b-a**2)

                #Calculating the Geometric integrals Mx & My
                Mx_calculated = (c_x/2)*(np.log((s_length[k]**2 + 2*a*s_length[k] + b)/(b))) + ((d_x-a*c_x)/(e))*(np.arctan2((s_length[k]+a),(e)) - np.arctan2((a),(e)))
                My_calculated = (c_y/2)*(np.log((s_length[k]**2 + 2*a*s_length[k] + b)/(b))) + ((d_y-a*c_y)/(e))*(np.arctan2((s_length[k]+a),(e)) - np.arctan2((a),(e)))
                
                #Checking for NANINF , Division by zero , or Imaginary
                if (np.iscomplex(Mx_calculated) or np.isnan(Mx_calculated) or np.isinf(Mx_calculated)):        
                    Mx_calculated = 0                                                           
                if (np.iscomplex(My_calculated) or np.isnan(My_calculated) or np.isinf(My_calculated)):         
                    My_calculated = 0                   

                Mx[i,j,k] =  Mx_calculated
                My[i,j,k] =  My_calculated

            
    return x_grid , y_grid , Mx , My

#Calculating the Vector fields 
def sp_calculate_vector_fields( x_grid , y_grid ,Mx , My , alpha , v_inf , x_bound , y_bound , panel_strengths):
    
    # initialize vector fields and create the path for the next segment
    Vx = np.zeros_like(x_grid)
    Vy = np.zeros_like(y_grid)
    coords = np.column_stack((x_bound , y_bound))
    path_c = path.Path(coords)

    #Determine if the point is in the geometry
    xlist = x_grid.shape
    ylist = x_grid.shape
    for i in range(ylist[0]):
        for j in range(xlist[1]):
            if path_c.contains_point([x_grid[i,j],y_grid[i,j]]):
                Vx[i,j] = 0
                Vy[i,j] = 0
            else:
                Vx[i,j] = v_inf*np.cos(alpha) + np.sum((panel_strengths*Mx[i,j,:])/(2*np.pi))
                Vy[i,j] = v_inf*np.sin(alpha) + np.sum((panel_strengths*My[i,j,:])/(2*np.pi))
                
    return Vx,Vy

#Calculating pressure distributions
def sp_calculate_pressure_distribution(J,v_inf,panel_strengths,beta, x_cont):
    vt = np.zeros_like(panel_strengths)
    C_p = np.zeros_like(panel_strengths)
    C_p_ref = np.zeros_like(panel_strengths)
    for i in range(len(vt)):
        vt[i] = v_inf*np.sin(beta[i]) + sum(J[i,:]*panel_strengths[i]/(2*np.pi))
        C_p[i] = 1 + ((vt[i])/v_inf)**2

    C_p_ref = 1 - 4*np.sin(x_cont)**2

    return C_p , C_p_ref
  
