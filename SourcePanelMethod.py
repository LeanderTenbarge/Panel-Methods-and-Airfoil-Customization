#Leander Tenbarge
# Source Panel formulation for non-lifting cylinder
import SourcePanelFunctions
import numpy as np
from matplotlib import pyplot

#Input Parameters
alpha = 0
v_inf = 1
x_dim = np.array([-5,5])
y_dim = np.array([-3,3])
num_points = 75
numB = 50
tO   = (360/(numB-1))/2
theta = np.linspace(0,360,numB)                                                 # Create angles for computing boundary point locations [deg]
theta = theta + tO                                                              # Add panel angle offset [deg]
theta = theta*(np.pi/180)                                                       # Convert from degrees to radians [rad]

# Boundary points
x = np.cos(theta)                                                              # Compute boundary point X-coordinate [radius of 1]
y = np.sin(theta)  

#calling the functions:
[x_bound, y_bound] = SourcePanelFunctions.check_direction(x,y)
[x_cont, y_cont , s_length , phi] = SourcePanelFunctions.process_geometry(x_bound, y_bound)
[I,J] = SourcePanelFunctions.sp_geometric_integrals(x_cont,y_cont,x_bound,y_bound,phi,s_length)
[beta, panel_strengths] = SourcePanelFunctions.find_source_panel_strengths(I , phi , alpha , v_inf)
[x_grid , y_grid , Mx , My ] = SourcePanelFunctions.sp_calc_point_geometric(num_points,x_dim,y_dim,x_bound,y_bound , phi ,s_length)
[Vx,Vy] = SourcePanelFunctions.sp_calculate_vector_fields( x_grid , y_grid ,Mx , My , alpha , v_inf , x_bound , y_bound , panel_strengths)
[C_p , C_p_ref]=SourcePanelFunctions.sp_calculate_pressure_distribution(J,v_inf,panel_strengths,beta,x_cont)

#Plotting the Results
pyplot.plot(x_bound, y_bound)
y_stream = np.linspace(y_dim[0] , y_dim[1] , 75)
x_stream = np.ones_like(y_stream)*x_dim[0]
start_stream_points= np.vstack((x_stream , y_stream)).T
pyplot.streamplot(x_grid,y_grid,Vx,Vy , density = 5 , color = "blue" , linewidth= .5 , arrowstyle="-" )
pyplot.plot(x_cont,C_p, label = 'Calculated')
#pyplot.quiver(x_grid,y_grid,Vx,Vy)
pyplot.show()

