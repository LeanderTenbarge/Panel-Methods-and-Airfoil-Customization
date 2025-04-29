#Leander Tenbarge Vortex Panel Method
import SourcePanelFunctions
import VortexPanelFunctions
import AirfoilTools
import numpy as np
from matplotlib import pyplot


# Input Parameters
alpha = .17
v_inf = 10
num_points = 100
x_dim = np.array([-2,2])
y_dim = np.array([-1.5,1.5])


#Getting the airfoil coordinates from the .txt file
[x_coords , y_coords] = AirfoilTools.get_coordinates('S1020.txt')

#Processing the geometry
[x_bound, y_bound] = SourcePanelFunctions.check_direction(x_coords,y_coords)
[x_cont, y_cont , s_length , phi] = SourcePanelFunctions.process_geometry(x_bound, y_bound)

# Calculating the Geometric integrals
[K,L] = VortexPanelFunctions.vp_geometric_integrals(x_cont,y_cont,x_bound,y_bound,phi,s_length)
[beta,panel_strengths] = VortexPanelFunctions.find_vortex_panel_strengths(K,phi,alpha,v_inf)
[x_grid,y_grid,Nx,Ny] = VortexPanelFunctions.vp_calc_point_geometric(num_points, x_dim, y_dim, x_bound , y_bound,phi,s_length)
[Vx,Vy] = VortexPanelFunctions.vp_calculate_vector_fields(x_grid, y_grid, Nx, Ny, alpha, v_inf, x_bound, y_bound, panel_strengths)
[Vt,Cp] = VortexPanelFunctions.vp_pressure_distribution(beta,panel_strengths,L,v_inf)
VortexPanelFunctions.find_parameters(Vt , s_length , v_inf)


#Plotting

pyplot.plot(x_bound, y_bound)
#pyplot.plot(x_cont,Cp,"o")
y_stream = np.linspace(y_dim[0] , y_dim[1] , 75)
x_stream = np.ones_like(y_stream)*x_dim[0]
start_stream_points= np.vstack((x_stream , y_stream)).T
pyplot.streamplot(x_grid,y_grid,Vx,Vy , density = 10 , color = "blue" , linewidth= .5 , arrowstyle="-" )
#pyplot.quiver(x_grid,y_grid,Vx,Vy)
pyplot.show()


