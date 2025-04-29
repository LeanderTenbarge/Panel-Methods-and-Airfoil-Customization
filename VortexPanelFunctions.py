#Leander Tenbarge 
#Vortex panel functions, the functions for geometry processing and manipulation are in the source panel functions script.
# Most of the Functions are extremely similar to the source panel aside from key differences (Kutta Condition, Etc.....)
import numpy as np
from matplotlib import path
import math


def vp_geometric_integrals(x_cont,y_cont,x_bound,y_bound,phi,s_length):
    #Initalizing the Geometric integral arrays
    K = np.zeros([len(x_cont),len(x_cont)])
    L = np.zeros([len(x_cont),len(x_cont)])

    #looping over the matricies
    for i in range(len(x_cont)):
        for j in range(len(x_cont)):
            if i != j:

                #Solving for the coefficients 
                a = - (x_cont[i] - x_bound[j])*np.cos(phi[j]) - (y_cont[i]-y_bound[j])*np.sin(phi[j])
                b = (x_cont[i]-x_bound[j])**2+(y_cont[i]-y_bound[j])**2
                c_l = np.sin(phi[j]-phi[i])
                d_l = (x_cont[i] - x_bound[j])* np.sin(phi[i]) - (y_cont[i]-y_bound[j])*np.cos(phi[i])
                c_k = - np.cos(phi[i]-phi[j])
                d_k = (x_cont[i]-x_bound[j])*np.cos(phi[i]) + (y_cont[i] - y_bound[j])*np.sin(phi[i])
                e = np.sqrt(b-a**2)

                #Finding the arrays of geometric integrals
                K[i,j] = (c_k/2)*(np.log((s_length[j]**2 + 2*a*s_length[j] + b)/(b))) + ((d_k-a*c_k)/(e))*(np.arctan2((s_length[j]+a),(e)) - np.arctan2((a),(e)))
                L[i,j] = (c_l/2)*(np.log((s_length[j]**2 + 2*a*s_length[j] + b)/(b))) + ((d_l-a*c_l)/(e))*(np.arctan2((s_length[j]+a),(e)) - np.arctan2((a),(e)))

                # Check if the value is NANinf or divide by zero to mitigate chance of numerical error
                if (np.iscomplex(K[i,j]) or np.isnan(K[i,j]) or np.isinf(K[i,j])):      
                    K[i,j] = 0                                                          
                if (np.iscomplex(L[i,j]) or np.isnan(L[i,j]) or np.isinf(L[i,j])):      
                    L[i,j] = 0             
   

    return K , L

def find_vortex_panel_strengths(K,phi,alpha,v_inf):

    #Creating the A Matrix and enforcing the kutta condition
    K_calc = -1*K
    K_calc[-1,:] = 0 
    K_calc[-1,0] = 1 
    K_calc[-1,-1] = 1

    #Constructing the beta(i) array 
    beta = np.zeros(len(K))
    for i in range(len(phi)):
        beta[i] = phi[i] + (np.pi/2) - alpha

    #Contructing the B matrix in x = a^-1*b
    b_calc = np.zeros(len(K))
    for i in range(len(K)):
        b_calc[i] = -v_inf * 2 * np.pi * np.cos(beta[i])

    #Enforcing the Kutta condition
    b_calc[-1] = 0 

    #solving the equation
    panel_strengths = np.zeros(len(K_calc))
    panel_strengths = np.linalg.solve(K_calc,b_calc)

    return beta  , panel_strengths

def vp_calc_point_geometric(N_vectors,x_dim,y_dim,x_bound,y_bound,phi ,s_length):

    #Initializing the arrays
    N_panels = int(len(x_bound) - 1 )
    x = np.linspace( x_dim[0] , x_dim[1] , int(((x_dim[1]-x_dim[0])/(y_dim[1]-y_dim[0])) * N_vectors))
    y = np.linspace( y_dim[0] , y_dim[1] , N_vectors)
    z = np.linspace(0 , N_panels-1 , N_panels)
    x_grid , y_grid = np.meshgrid(x,y,indexing='xy')
    Nx = np.zeros((len(y),len(x),len(z)))
    Ny = np.zeros((len(y),len(x),len(z)))
   

    # Creating the geomtric integral 
    for i in range(len(y)): # X points
        for j in range(len(x)): # Y points
            for k in range(len(z)):  # Iterating over every single panel for each point:
            
                #Calculating the Coeffcients 
                a = - (x_grid[i,j] - x_bound[k])*np.cos(phi[k]) - (y_grid[i,j] - y_bound[k])*np.sin(phi[k])
                b = (x_grid[i,j] - x_bound[k])**2 + (y_grid[i,j] - y_bound[k])**2
                c_x = np.sin(phi[k])
                c_y = -np.cos(phi[k])
                d_y = (x_grid[i,j] - x_bound[k])
                d_x = -(y_grid[i,j] - y_bound[k])
                e = np.sqrt(b-a**2)

                #Calculating the Geometric integrals Mx & My
                Nx_calculated = (c_x/2)*(np.log((s_length[k]**2 + 2*a*s_length[k] + b)/(b))) + ((d_x-a*c_x)/(e))*(np.arctan2((s_length[k]+a),(e)) - np.arctan2((a),(e)))
                Ny_calculated = (c_y/2)*(np.log((s_length[k]**2 + 2*a*s_length[k] + b)/(b))) + ((d_y-a*c_y)/(e))*(np.arctan2((s_length[k]+a),(e)) - np.arctan2((a),(e)))
                
                #Checking for NANINF , Division by zero , or Imaginary
                if (np.iscomplex(Nx_calculated) or np.isnan(Nx_calculated) or np.isinf(Nx_calculated)):        
                    Nx_calculated = 0                                                           
                if (np.iscomplex(Ny_calculated) or np.isnan(Ny_calculated) or np.isinf(Ny_calculated)):         
                    Ny_calculated = 0                   

                Nx[i,j,k] =  Nx_calculated
                Ny[i,j,k] =  Ny_calculated

            
    return x_grid , y_grid , Nx , Ny

def vp_calculate_vector_fields(x_grid, y_grid, Nx , Ny, alpha, v_inf, x_bound , y_bound , panel_strengths):
    # initialize vector fields and create the path for the next segment
    Vx = np.zeros_like(x_grid)
    Vy = np.zeros_like(y_grid)
    coords = np.column_stack((x_bound , y_bound))
    path_c = path.Path(coords)

    #Determine if the point is in the geometry
    xlist = x_grid.shape
    ylist = x_grid.shape

    #Determine if the point is in the geometry and assigning a velocity if not in the shape. 
    xlist = x_grid.shape
    ylist = x_grid.shape
    for i in range(ylist[0]):
        for j in range(xlist[1]):
            if path_c.contains_point([x_grid[i,j],y_grid[i,j]]):
                Vx[i,j] = 0
                Vy[i,j] = 0
            else:
                Vx[i,j] = v_inf*np.cos(alpha) + np.sum((-panel_strengths*Nx[i,j,:])/(2*np.pi))
                Vy[i,j] = v_inf*np.sin(alpha) + np.sum((-panel_strengths*Ny[i,j,:])/(2*np.pi)) 
                
    return Vx,Vy

def vp_pressure_distribution(beta,panel_strengths,L,v_inf):
    Vt = np.zeros_like(panel_strengths)
    Cp = np.zeros_like(Vt)
    for i in range(len(Vt)):
        Vt[i] = v_inf*np.sin(beta[i]) + panel_strengths[i]/2 + sum(((-panel_strengths[i]/(2*np.pi)))*L[i,:])
        Cp[i] = 1 - (Vt[i]/v_inf)**2

    return Vt , Cp

def find_parameters(Vt , s_length, v_inf):
    circ = sum(Vt[:]*s_length[:])
    lift = (1.225)*(v_inf)*circ
    Cl = lift/(.5*1.225*v_inf**2)
    print(lift,Cl)