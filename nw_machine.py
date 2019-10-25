import networkx as nx
import matplotlib.pyplot as plt
import cmath
import sys
import numpy as np
from numpy import pi, sin, exp, mean
import matplotlib.pyplot as plt
from networkx.algorithms import community
#from networkx.algorithms.community import LFR_benchmark_graph
import os

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

def SW(N,k,p,wseed):
	if wseed == 1:
		G = nx.watts_strogatz_graph(N, k, p, seed=123)
	else:
		G = nx.watts_strogatz_graph(N, k, p)
	return G

def Comm(N, wseed):

	if wseed == 1:
		np.random.seed(50)

	N = 1000
	tau1 = 2
	tau2 = 1.1
	mu = 0.05
	max_deg = 0.1*N
	k = 8
	minc = 30
	maxc = 0.1*N

	G = LFR_benchmark_graph(n = N, tau1 = tau1, tau2 = tau2, mu = mu, min_degree = k, max_degree = k, min_community=minc, max_community = maxc, seed = 10)
                        
	return G
	
def BA(N,m,wseed):
	if wseed == 1:
		G = nx.barabasi_albert_graph(N, m, seed=123)
	else:
		G = nx.barabasi_albert_graph(N, m)
	return G

def ER(N, avg_deg, wseed):
	if wseed == 1:
		G = nx.erdos_renyi_graph(N, avg_deg/N, seed=123, directed=False)
	else:
		G = nx.erdos_renyi_graph(N, avg_deg/N, directed=False)	
	return G
	
def RK4(w, K, phi, A, klist):	
	#algorithm function of RungeKutta of 4th order, input the initial vector phi, outputs the new vector phi, h the time step. 
	#phi, k1 to k4 are all N-Vector	
	k1 = h * ode(w, K, phi, A)
	k2 = h * ode(w, K, phi + k1 / 2e0, A)
	k3 = h * ode(w, K, phi + k2 / 2e0, A)
	k4 = h * ode(w, K, phi + k3, A)
	phi2 = phi + (k1 + 2e0 * k2 + 2e0 * k3 + k4) / 6e0
	return phi2

def sRK4(f, g, w, eps, var, A, eta, psi, sigma):
	#strong order 1.5 with Stratonovich interpretation
	odeK1 = f(sigma, w, eps, var, A)
	K2 = var + TT * h * odeK1 + TT * eta * g(var, sigma)
	odeK2 = f(sigma, w, eps, K2, A)
	K3 = var + TH * h * odeK1 - OT * h * odeK2 + OH * eta * g(var, sigma) + OS * eta * g(K2, sigma) - TT * psi * g(var, sigma)
	K4 = var + SS * h * odeK1 - OH * eta * g(var, sigma) + OH * eta * g(K3, sigma) + OS * psi * g(var, sigma) + OH * psi * g(K2, sigma)
	var = var + h * (OF * odeK1 + TF * odeK2 - TF * f(sigma, w, eps, K3, A) + TF * f(sigma, w, eps, K4, A)) + eta * (- OH * g(var, sigma) + TH * g(K2, sigma) - TF * g(K3, sigma) + TF * g(K4, sigma)) + psi * (TH * g(var, sigma) - TH * g(K2, sigma))
	return var

def sRK42(f, w, eps, var, A, eta, psi, sigma):
	#strong order 1.5 with Stratonovich interpretation for individual noise
	odeK1 = f(sigma, w, eps, var, A)
	K2 = var + TT * h * odeK1 + TT * eta * sigma
	odeK2 = f(sigma, w, eps, K2, A)
	K3 = var + TH * h * odeK1 - OT * h * odeK2 + OH * eta * sigma + OS * eta * sigma - TT * psi * sigma
	K4 = var + SS * h * odeK1 - OH * eta * sigma + OH * eta * sigma + OS * psi * sigma + OH * psi * sigma
	var = var + h * (OF * odeK1 + TF * odeK2 - TF * f(sigma, w, eps, K3, A) + TF * f(sigma, w, eps, K4, A)) + eta * (- OH * sigma + TH * sigma - TF * sigma + TF * sigma) + psi * (TH * sigma - TH * sigma)
	return var
			
def g1(var, sigma):
	return sigma * sin(var)
	
def ode1(sigma, w, eps, phi, A):
	return w + K *( A.dot(exp(1j * phi)) * exp(-1j * phi) ).imag + sigma ** 2e0 * sin(2e0 * phi)

def ode2(sigma, w, eps, phi, A):
	return w + K *( A.dot(exp(1j * phi)) * exp(-1j * phi) ).imag + sin(phi)

def ode3(sigma, w, eps, phi, A):
	return w + sin(phi) * K *( A.dot(exp(1j * phi)) * exp(-1j * alpha) ).imag + sigma ** 2e0 * sin(2e0 * phi)

def ode4(sigma, w, eps, phi, A):
	return w + K *( A.dot(exp(1j * phi)) * exp(-1j * phi) ).imag
				
def C_k_local(N, var, klist):
	sumN = np.dot(klist, np.exp(1j * var))
	R, Theta = cmath.polar(sumN)
	return R/np.sum(klist), Theta

def C_k(N, var):
	sumN = np.sum(np.exp(1j * var))
	R, Theta = cmath.polar(sumN)
	return R/N, Theta
		
def integrate(A, klist, K, N, phi, w, sigma, eqn, datfile):
	if eqn == '1':
		for j in range(intN):
			u = np.random.normal(0, 1, N)
			v = np.random.normal(0, 1, N)
			eta = u * sqrth
			psi = sqrth * (u/2e0 + v/c)	
		
			phi = sRK4(ode1, g1, w, K, phi, A, eta, psi, sigma) # multiplicative nosie

	if eqn == '4_indnoise':
		for j in range(intN):
			u = np.random.normal(0, 1)
			v = np.random.normal(0, 1)
			eta = u * sqrth
			psi = sqrth * (u/2e0 + v/c)	
		
			phi = sRK42(ode4, w, K, phi, A, eta, psi, sigma) # multiplicative nosie
						
	R_local = C_k_local(N, phi, klist)[0]
	R_global = C_k(N, phi)[0]
	
	f = open(datfile, 'a')
	f.write(str(R_local) + '\t' + str(R_global) +'\n')
	f.close()

def createNWs():
	for i, NW in enumerate(NW_type):
		for k in range(numofNWs):
			if NW == 'SF':
				G = Scale_free(N, 2.3, 0)
			elif NW == 'BA':
				G = BA(N, 5, 0) # scale free Barabási-Albert preferential attachment network #avg deg = 50 if BA_parameter = 25
			elif NW == 'ER':
				G = ER(N, 10, 0) # 3rd arg = with seed (1) or not (0)	
			elif NW == 'WS001':
				G = SW(N, 10, 0.01, 0) 
			elif NW == 'WS01':
				G = SW(N, 10, 0.1, 0) 
			elif NW == 'Comm':
				G = Comm(N, 0) 
		#			elif NW == 'Asso':
			
			while nx.is_connected(G) == 0:
				num_of_disconn_graphs = nx.number_connected_components(G)
				nodes_tobeconn = []
				sortedsubgraph = sorted(nx.connected_components(G), key = len, reverse=True)
				for i in range(num_of_disconn_graphs):
					disconn_subgr = list(sortedsubgraph[i])
					nodeid_insubgr = np.random.randint(len(disconn_subgr))
					nodes_tobeconn.append(disconn_subgr[nodeid_insubgr])
				for j in range(len(nodes_tobeconn)-1):
					G.add_edge(nodes_tobeconn[j], nodes_tobeconn[j+1])

			fname = NW + str(k) + 'N=' + str(N) + '.dat'
			nx.write_adjlist(G, os.path.join('NWs', fname))
	sys.exit()
	
if __name__ == "__main__":
	N = 500
	hlist = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2]
#	NW_type = ['BA', 'ER', 'WS001', 'WS01', 'COMM']
	NW_type = ['BA', 'ER']
	eqn = '1'
	
	c = 2e0 * np.sqrt(3e0)
	
	sigma_list = np.arange(0, 1e0, 1e-1)	
	
	TT = 2e0 / 3e0
	OT = 1e0 / 3e0
	TH = 3e0 / 2e0
	OH = 1e0 / 2e0
	OS = 1e0 / 6e0
	SS = 7e0 / 6e0
	TF = 3e0 / 4e0
	OF = 1e0 / 4e0
	
	phi = 2e0 * pi * np.random.rand(N)
		
	Nrun = 8

	Klist = [2e-1, 3e-1, 4e-1, 5e-1]
	
#	w = 10e0 * np.random.rand(N)
	w = np.random.standard_cauchy(1000000)
	w = w[(w>-10) & (w<10)]
#	print(len(w))

	interval = len(w)//500 
	w = w[:-interval][::interval]
	
#	print(len(w))
	numofNWs = 50
	
#	createNWs()
	
	processornum = int(sys.argv[1])-1
	noiserealization = 1
	run = processornum%noiserealization

	#	run = 1
	for K in Klist: 
		for i, NW in enumerate(NW_type):
			h = hlist[i]
			intN = int(5e0 / h)
			sqrth = np.sqrt(h)
	#		for sigma in sigma_list:
			sigma = sigma_list[processornum//noiserealization]
			datfile = 'R_noise' + str(run) + '_sigma_=' + str(sigma) + '_sRK4_' + NW +'_SDEtype' + eqn + 'cauchy_K=' + str(K) + '.dat'
			f = open(datfile, 'w')
			f.close()
			for k in range(numofNWs): #number of NWs

				fname = NW + str(k) + 'N=' + str(N) + '.dat'
				G = nx.read_adjlist(os.path.join('NWs', fname))
	
				A = nx.to_numpy_matrix(G)
				A = np.array(A.tolist())
				klist = np.sum(A, axis = 0)

				integrate(A, klist, K, N, phi, w, sigma, eqn, datfile)
		
