import numpy as np
import matplotlib.pyplot as pyp
import pylab as pyl
import math
import timeit

begin = timeit.default_timer()	# tic

pi = 3.14						# pi

# simulation variable
n_reach = 101					# number of segments
n_itr = 5001					# number of iterations
stat_val = "y"

# pipe properties
L = 1000						# length of pipe (m)
D = 0.634						# diameter of pipe (m)
E = 2.4 * math.pow(10,9)		# modulus of elasticity of PVC (N/m^2)
e_t = 0.0132					# wall thickness (m)
f = 0.02						# Darcy-Weisbach friction factor
ke = 0							# coefficient of entrance head loss of reservoir
theta = math.atan(0)			# pipe's declination (rad), arg: altitude/distance
A = 0.25 * pi * math.pow(D,2)	# pipe's cross section area (m^2)
R = f / (2 * D * A)				# friction constant for char equation

# water properties
K = 2.15 * math.pow(10,9)		# bulk modulus (N/m^2)
rho = 998.2						# mass density at 20 deg C
a = 250							# wave velocity (m/s)

# general properties
g = 9.81						# gravity acceleration (m/s^2)
Q_ss = 0.1						# steady state discharge rate (m/s^3)
H_res = 3						# steady state pressure head (m)
tau = 1							# state of valve; 1: fully open

delta_x = L / (n_reach - 1)					# length of reach (m)
delta_t = L / ((n_reach - 1) * a)			# sampling period (s)
t_max = math.ceil((n_itr - 1) * delta_t)	# maximum simulation time (s)

Q = np.zeros((n_itr,n_reach))	# define matrix for discharge rate
H = np.zeros((n_itr,n_reach))	# define matrix for pressure head

Ca = g * A / a						# coefficient a
Cp = np.zeros((n_itr,n_reach-1))	# initialisation of positive char. eq matrix
Cn = np.zeros((n_itr,n_reach-1))	# initialisation of negative char. eq matrix
Cv = np.zeros((n_itr,1))			# initialisation of valve state matrix

t = np.arange(0,t_max + delta_t,delta_t)	# range of t
n_counter = 1		# counter

for i in range(n_reach):
	Q[0][i] = Q_ss				# initialisation of discharge rate along pipe
	H[0][i] = H_res - (i * f * delta_x * math.pow(Q_ss,2) / (2 * g * D * A)) 	# initialisation of pressure head along pipe
	
while n_counter < n_itr:
	
	# calculation of characteristic equation parameter
	for i in range(n_reach-1):
		Cp[n_counter][i] = Q[n_counter-1][i] + (Ca * H[n_counter-1][i]) - (R * delta_t * Q[n_counter-1][i] * abs(Q[n_counter-1][i]))
		Cn[n_counter][i] = Q[n_counter-1][i+1] - (Ca * H[n_counter-1][i+1]) - (R * delta_t * Q[n_counter-1][i+1] * abs(Q[n_counter-1][i+1]))
	
	# at upstream boundary
	H[n_counter][0] = H_res
	Q[n_counter][0] = Cn[n_counter][0] + (Ca * H[n_counter][0])
	
	# along pipe
	for i in range(1,n_reach-1):
		Q[n_counter][i] = 0.5 * (Cp[n_counter][i-1] + Cn[n_counter][i])
		H[n_counter][i] = (Cp[n_counter][i-1] - Q[n_counter][i]) / Ca
		
	# at downstream boundary
	if stat_val == 'y' and n_counter == 1:
		tau = 0
	
	Cv[n_counter][0] = tau * math.pow(Q[0][n_reach-1],2) / (Ca * H[n_counter-1][n_reach-1])
	Q[n_counter][n_reach-1] = 0.5 * (-Cv[n_counter][0] + math.sqrt(math.pow(Cv[n_counter][0],2) + (4 *  Cv[n_counter][0] * Cp[n_counter][n_reach-2])))
	H[n_counter][n_reach-1] = (Cp[n_counter][n_reach-2] - Q[n_counter][n_reach-1]) / Ca
	
	n_counter += 1

# plotting pressure head and discharge rate
pyp.plot(t,H[:,n_reach-1])
pyp.xlim([0,t_max])
pyp.xlabel('Time (s)')
pyp.ylabel('Pressure Head (m)')
pyp.title('Pressure Head at Downstream Valve')
pyp.grid(True)

stop = timeit.default_timer()		# toc
print(stop - begin)

pyp.show()
#just a comment