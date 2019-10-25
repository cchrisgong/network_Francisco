import networkx as nx
import matplotlib.pyplot as plt
import cmath
import sys
import numpy as np
from numpy import pi, sin, cos, exp
import matplotlib.pyplot as plt
from time import sleep

def wrap2(var):
	#returns angle (vector) var to interval [-pi,pi]
	return np.mod(var-pi,2e0*pi)-pi

def Scale_free(N, gamma, wseed):
	s = nx.utils.powerlaw_sequence(N, gamma) #100 nodes, power-law exponent 2.5
	if wseed == 1:
		G = nx.expected_degree_graph(s, seed=123, selfloops=False)
	else:
		G = nx.expected_degree_graph(s, selfloops=False)
	return G
	
def ER(N, avg_deg, wseed):
	if wseed == 1:
		G = nx.erdos_renyi_graph(N, avg_deg/N, seed=123, directed=False)
	else:
		G = nx.erdos_renyi_graph(N, avg_deg/N, directed=False)	
	return G
	
def RK4(omega, K, phi, A):	
	#algorithm function of RungeKutta of 4th order, input the initial vector phi, outputs the new vector phi, h the time step. 
	#phi, k1 to k4 are all N-Vector	
	k1 = h * ode(omega, K, phi, A)
	k2 = h * ode(omega, K, phi + k1 / 2e0, A)
	k3 = h * ode(omega, K, phi + k2 / 2e0, A)
	k4 = h * ode(omega, K, phi + k3, A)
	dphi = phi + (k1 + 2e0 * k2 + 2e0 * k3 + k4) / 6e0
	return dphi

def C_k(N, var, kls, order):
	sumN = np.dot(kls, np.exp(1j * order * var))
	R, Theta = cmath.polar(sumN)
	return R/np.sum(kls), Theta

def ode(omega, eps, phi, A):
	dphi = omega + eps *( A.dot(exp(1j * phi)) * exp(-1j * phi) ).imag
	return dphi

def g(var):
#	return sigma * sin(var) * omega
#	return sigma * sin(var)
	return sigma

def ode_strat(N, omega, eps, phi, A):
	return ode(omega, eps, phi, A) - strat_shift(phi)
	
def strat_shift(var):
	return sigma ** 2e0 * sin(2e0 * var) # 1/2 * (sigma ** 2e0) * sin(var)cos(var)
			
def integrate(G,K,N):
	while nx.is_connected(G) == 0:
		num_of_disconn_graphs = nx.number_connected_components(G)
		nodes_tobeconn = []
		for i in range(num_of_disconn_graphs):
			disconn_subgr = list(sorted(nx.connected_components(G), key = len, reverse=True)[i])
			nodeid_insubgr = int(np.random.rand() * int(len(disconn_subgr)))
			nodes_tobeconn.append(nodeid_insubgr)
		for j in range(len(nodes_tobeconn)-1):
			G.add_edge(nodes_tobeconn[j],nodes_tobeconn[j+1])
	
	A = nx.to_numpy_matrix(G)
	A = np.array(A.tolist())
#	A = np.ones((N,N))
#	pos = nx.spring_layout(G)
#	nx.draw_networkx(G, pos)
#	plt.show()
#	sys.exit()
	kls = np.sum(A, axis = 0)
	
	phi = phiin
	phiB = phi

#	B: deterministic
	for j in range(intN):
		u = np.random.normal(0, 1)
		v = np.random.normal(0, 1)
		eta = u * sqrth
		psi = sqrth * (u/2e0 + v/c)	
		
		phi2 = Euler_Heun(phi, h, sqrth, psdelta, xi)
#		phi2 = sRK4(ode_strat, omega, K, phi, A, eta, psi) # multiplicative nosie
		phi2B = RK4(omega, K, phiB, A)
		phiB = wrap2(phi2B)
		phi = wrap2(phi2)
		if j % 1000 == 0:
			time =j * h
			print(time, C_k(N, phi, kls, 1e0)[0])
			print(time, C_k(N, phiB, kls, 1e0)[0])			
#			print(phi[:10])
#			plotkura(phi, time, 0)
#			plotkura(phiB, time, 1)			
	R = C_k(N, phi, kls, 1e0)[0]
	R2 = C_k(N, phi, kls, 2e0)[0]
	RB = C_k(N, phiB, kls, 1e0)[0]
	R2B = C_k(N, phiB, kls, 2e0)[0]
	print(R, R2, RB, R2B)	
	return R, R2, RB, R2B
#	plt.scatter(omega, wrap2(phi), s = 1)
#	plt.ylabel('phi')
#	plt.xlabel('omega')
#	plt.show()
#	sys.exit()

def plotkura(phi, time, tag):
	# plotting the phase on a unit circle
	fig = plt.figure(2)
	fig.clf()
	ax = fig.add_subplot(111, aspect='equal', xlim=(-1.1,1.1), ylim=(-1.1,1.1))
	
	circ = plt.Circle((0, 0), radius=1, edgecolor='b', facecolor='None')
	ax.add_patch(circ)
	for j in range(len(phi)):
		if tag == 0:
			color0 = 'b'
		else:
			color0 = 'r'
		plt.scatter(cos(phi)[j], sin(phi)[j], color = color0)
#		ax.annotate(str(j), xy=(cos(phi)[j], sin(phi)[j]), textcoords='offset points')
	plt.grid()
	plt.show(block=False)
	sleep(2)
#	fig.savefig(os.path.join(os.getcwd(), 'case'+str(casenum), str(time)+"kura.png"), dpi=500)
	plt.close()
	
def plot():
	Rllist = []
	for run in range(Nrun):
#		datfile = 'K_vs_R_' + str(run) + '_sigma_=' + str(sigma) + '_alpha=' + str(psdelta) + '_sde_ER.dat'
		datfile = 'K_vs_R_' + str(run) + '_sigma_=' + str(sigma) + '_alpha=' + str(psdelta) + '_ode_ER.dat'

		K, R, R2 = np.loadtxt(datfile, unpack =1)
		Rllist.append(R)
	print(len(np.mean(Rllist, axis=0)),len(Klst))
	fig = plt.figure()
	ax = fig.add_subplot(111)

	ax.set_xlabel("K",fontsize=20)
	ax.set_ylabel("R", fontsize=20)	
	ax.plot(Klst, np.mean(Rllist, axis=0))
	
	fig.savefig('K_vs_R.jpg', dpi=500)
	plt.close()	
	
if __name__ == "__main__":

############### nonidentical (uniform frequency distn) with strat shift, lead quickly to 2 clusters? sRK4, repulsive
	Nrun = 7
	run = int(sys.argv[1]) - 1
	N = 100
	Klst = np.arange(1e-2, 1.1e0, 5e-2)
#	K = 1.1e0
		
	h = 1e-3
	sqrth = np.sqrt(h)
	intN = int(200e0 / h)
	c = 2e0 * np.sqrt(3e0)

#	psdelta = 0.3e0*2e0*pi
	psdelta = 0
#	psdelta = 0.5e0 * 2e0 * pi
	sigma = 1e-1
	
	TT = 2e0 / 3e0 
	OT = 1e0 / 3e0
	TH = 3e0 / 2e0 
	OH = 1e0 / 2e0
	OS = 1e0 / 6e0
	SS = 7e0 / 6e0
	TF = 3e0 / 4e0
	OF = 1e0 / 4e0
	
#	Graph = Scale_free(N, 2.5e0, 1)
	Graph = ER(N, 5e0, 0) # 3rd arg = with seed (1) or not (0)

	R_list = []
	R2_list = []
	RB_list = []
	R2B_list = []	
	for K in Klst:
		omega = 2e0 * np.random.rand(N) - 1e0
	#		omega = 0
		phiin = 2e0 * pi * np.random.rand(N)
		R_sto, R2_sto, R_det, R2_det = integrate(Graph, K, N)
		R_list.append(R_sto)
		R2_list.append(R2_sto)
		RB_list.append(R_det)
		R2B_list.append(R2_det)
		
	DAT = np.column_stack((Klst, R_list, R2_list))
	datfile1 = 'K_vs_R_' + str(run) + '_sigma_=' + str(sigma) + '_alpha=' + str(psdelta) + '_sde_ER.dat'
	np.savetxt(datfile1, DAT, delimiter=" ")
	
	DAT = np.column_stack((Klst, RB_list, R2B_list))
	datfile2 = 'K_vs_R_' + str(run) + '_sigma_=' + str(sigma) + '_alpha=' + str(psdelta) + '_ode_ER.dat'
	np.savetxt(datfile2, DAT, delimiter=" ")
	
#	plot()
