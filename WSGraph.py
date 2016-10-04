###################################################################
#CS6824 Homework1
#Zhen Guo
#10/2/16
#functions for Watts-Strogatz graphs
###################################################################

import numpy as np
import networkx as nx
import pylab as plt
import math
import random
import time

#WS-1 function to generate Watts Strogatz graph
def generate_watts_strogatz_graph(n, k, p):
    """randomised ring graphs"""
    ringG = ring_graph(n, k)
    graphSW = rewire(ringG, n, k, p)
    return graphSW 

def ring_graph(n, k):
    """create undirected, regular ring graph with nk/2 edges"""
    graph = np.array([x[:] for x in [[0]*n]*n])
    G = nx.from_numpy_matrix(graph)
    for i in range(n):
        for j in range(i + 1, i + 1 + math.ceil(k/2)):
            G.add_edge(*(i, j % n))
    return G

def rewire(ringG, n, k, p):
    """rewire every edge with probability p"""
    for i in range(math.ceil(k/2)):
        for j in range(n):
            flag = randflag(p)
            if flag == 1:  #decide to rewire this edge
                ncnodes = non_connected_nodes(ringG, j)
                lennodes = len(ncnodes)
                newnode = ncnodes[random.randint(0, lennodes - 1)]
                ringG.remove_edge(*(j, (j + i + 1)% n))
                ringG.add_edge(*(j, newnode))
    return ringG

def randflag(p):
    """determine if an edge will be rewired at the probability of p"""
    chance = random.random()
    if chance < p:
        return 1
    else:
        return 0

def non_connected_nodes(G, node):
    """generate all nodes not connected to x"""
    nodelist = G.nodes()
    neighbor = G.neighbors(node)
    for nei in neighbor:
        nodelist.remove(nei)
    nodelist.remove(node)
    return nodelist

#WS-2 average shortest path length and clustering coefficient
def average_shortest_path_length(G):
    """connected graph only"""
    len = nx.average_shortest_path_length(G)
    return len

def clustering_coefficient(G):
    """undirected graph only"""
    cluster = nx.clustering(G)
    coeff = sum(cluster.values()) / len(cluster)
    return coeff

#WS-3 p values
def p_range(k):
    p = []
    x = math.pow(10, -4.0 / k)
    for i in range(k, 0, -1):
        p.append(math.pow(x, i))
    return p

#WS-3 run test with different p values
def test_p_n(nlist, k, plist, f):
    runningtime = {}
    for n in nlist:
        time1 = time.time()
        f.write("n = %d, k = 20:\n" %n)
        for p in plist:
            G = generate_watts_strogatz_graph(n, k, p) 
            lp = average_shortest_path_length(G)
            cp = clustering_coefficient(G)
            f.write("p = %.3e, l = %.2f, c = %.2f \n" %(p, lp, cp))   
        time2 = time.time()
        runningtime[n] = time2 - time1
    for n in nlist:
        f.write("Time used for running 20 times n = %d is: %.4fs \n" %(n, runningtime[n]))    
    print("WS-3 p tests completed.")
    print("Time used: %.4fs" %sum(runningtime.values()))
    return 



#WS-4 generate graph and output
def SW_figure1_plot(n, k, prange):
    """compute 20 p values with 10 runs, plot the mean and std of l and c"""
    time0 = time.time()
    p0 = 0
    p1 = 1
    #regular lattice
    G = generate_watts_strogatz_graph(n, k, p0) 
    l0 = average_shortest_path_length(G)
    c0 = clustering_coefficient(G)

    #rewire the regular lattice
    lave = [] #l values for 20 p
    cave = [] #c values for 20 p
    for p in prange:
        llist = []
        clist = []
        for i in range(10):     #run each p value 10 times
            G = generate_watts_strogatz_graph(n, k, p) 
            l = average_shortest_path_length(G)
            c = clustering_coefficient(G)
            llist.append(l / l0)
            clist.append(c / c0)
        npllist = np.array(llist)
        npclist = np.array(clist)
        lave.append([npllist.mean(), npllist.std()])
        cave.append([npclist.mean(), npclist.std()])

    #random graph
    l1list = []
    c1list = []
    for i in range(10):
        Grand = generate_watts_strogatz_graph(n, k, p1)
        l1 = average_shortest_path_length(Grand)
        c1 = clustering_coefficient(Grand)
        l1list.append(l1 / l0)
        c1list.append(c1 / c0)
    npl1 = np.array(l1list)
    npc1 = np.array(c1list)
    l1ave = [npl1.mean(), npl1.std()]
    c1ave = [npc1.mean(), npc1.std()]

    #output results to file
    time1 = time.time()
    runningtime = time1 - time0
    path = 'results/WS4_n_{}.out'.format(n)
    f = open(path, 'w')
    f.write("p = 0, l0 = %.3f, c0 = %.3f \n" %(l0, c0))
    for p in prange:
        f.write("p = %.3e, l/l0 = %.3f %.3e, c/c0 = %.3f %.3e\n" \
            %(p, lave[prange.index(p)][0],lave[prange.index(p)][1], \
            cave[prange.index(p)][0], cave[prange.index(p)][1]))
    f.write("p = 1, l1/l0 = %.3f %.3e, c1/c0 = %.3f %.3e \n" \
        %(l1ave[0], l1ave[1], c1ave[0], c1ave[1]))
    f.write("Time used: %.4fs" %runningtime)
    f.close()

    #draw graph based on the result
    nplave = np.array(lave)
    npcave = np.array(cave)
    lave_mean = nplave[:,0]
    lave_std = nplave[:,1]
    cave_mean = npcave[:,0]
    cave_std = npcave[:,1]
    lave_mean = np.insert(lave_mean, 0, 1)
    lave_std = np.insert(lave_std, 0, 0)
    cave_mean = np.insert(cave_mean, 0, 1)
    cave_std = np.insert(cave_std, 0, 0)
    lave_mean = np.append(lave_mean, l1ave[0])
    lave_std = np.append(lave_std, l1ave[1])
    cave_mean = np.append(cave_mean, c1ave[0])
    cave_std = np.append(cave_std, c1ave[1])

    ax = plt.subplot()
    ax.set_xscale("log", nonposx='clip')
    npp = np.array(prange)
    npp = np.insert(npp, 0, 0.00008)
    npp = np.append(npp, 1)
    plt.ylim([0, 1.05])
    plt.xlim([0.00008, 1])
    plt.errorbar(npp, lave_mean, yerr = lave_std, linestyle = 'None', marker = 's', markersize = 4, label = 'l(p)/l(0)')
    plt.errorbar(npp, cave_mean, yerr = cave_std, linestyle = 'None', marker = 'o', markersize = 4, label = 'c(p)/c(0)')
    title = 'n = {}, k = {}'.format(n, k)
    ax.set_title(title)
    plt.legend()
    plt.xlabel('p')
    path = 'results/WS4_n_{}.png'.format(n)
    plt.savefig(path) # save as png
    plt.clf()
    print("WS4 n={}, k=20 plot finished.".format(n))
    print("Time used: %.4fs" %runningtime)
    return