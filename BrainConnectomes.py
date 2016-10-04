###################################################################
#CS6824 Homework1
#Zhen Guo
#10/2/16
#functions for Analasis of Brain Connectomes
#I borrowed parse mat codes from Jeff Law but I deleted node names.
###################################################################

import WSGraph
import scipy
from scipy.io import loadmat
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
from optparse import OptionParser
import math
import random
import time

# format of the matlab files
MLFORMAT = {
    "primate_cortex": 'CIJ',
    "macaque_cortex": 'CIJ',
    "cat_cortex": 'CIJctx',
    "macaque_functional": 'CIJ',
    "coactivation": 'Coactivation_matrix',
}
YEARS = {
    "primate_cortex": "1991",
    "macaque_cortex": "1993",
    "cat_cortex": "1995",
    "macaque_functional": "2007",
    "coactivation": "2013",
    "high_density": "2012",
}

#BC-1 parsing files
def parse_mat_files(mlfile, path):
    """read network and return a networkx object"""
    print("\n" + '-'*20)
    print("loading %s (%s)" %(mlfile, path))
    loaded_mlfile = loadmat(path)
    ml_matrix = loaded_mlfile[MLFORMAT[mlfile]]
    G = matlab_to_nx(ml_matrix, name=mlfile)
    if G == -1:
        # if the parser failed, skip this entry
        print("Quitting")
        sys.exit()

    out_file = 'results/' + path.replace('.mat','') + '.txt'
    print ("writing %s to: %s" %(mlfile, out_file))
    with open(out_file, 'w') as out:
        if mlfile == "coactivation": #coactivation is undirected
            out.write('#node1\tnode2\tedge_weight\n')
        else:
            out.write('#tail\thead\tedge_weight\n')
        out.write('\n'.join(['%s\t%s\t%s' % (u,v,str(G.edge[u][v]['weight'])) for u,v in sorted(G.edges())]))
    return G

def matlab_to_nx(ml_matrix, name=''):
    """ function to parse a scipy.io.loadmat loaded matlab file and store it as a networkx graph object
    *ml_matrix*: matlab (numpy) matrix
    *name*: name of the graph (optional)
    *returns*: a networkx graph object or -1 if (# rows != # of columns or # rows != # names)
    """
    print("# of rows and columns in ml_matrix: %d %d" % (len(ml_matrix), len(ml_matrix[0])))
    if len(ml_matrix) != len(ml_matrix[0]):
        print("Error: # of rows in ml_matrix != # columns")
        return -1
    if name == "coactivation": 
        G = nx.Graph(name=name) #coactivation graph is undirected
    else: 
        G = nx.DiGraph(name=name)
    for i in range(len(ml_matrix)):
        for j in range(len(ml_matrix[i])):
            if ml_matrix[i][j] > 0:
                u = i
                v = j
                G.add_edge(u, v, weight=ml_matrix[i][j])

    print (nx.info(G))
    return G

def parse_high_density_file(mlfile, path):
    """read high-density graph and return a networkx object"""
    print ("\n" + '-'*20)
    print ("loading %s (%s)" %(mlfile, path))
    f = open(path, 'r')
    text = f.readlines()
    ml_matrix = []   #G91*29
    for line in text[1:]:
        datas = line.split()
        ml_matrix.append(datas[2:5])  #SOURCE TARGET FLNe
    f.close()
    #construct G29*29
    np_matrix = np.array(ml_matrix)
    nodelist = np.unique(np_matrix[:,1]) #29 nodes
    deletenodes = []
    for i in range(len(np_matrix)):
        if not np_matrix[i][0] in nodelist:
            deletenodes.append(i)
    np_matrix = np.delete(np_matrix, deletenodes, 0) #delete 62 nodes
    matrixlength = len(np_matrix)
    print(matrixlength)
    edgesrepeats = np.ones(matrixlength) #record edge repeat numbers
    G = nx.DiGraph(name=mlfile)
    for i in range(matrixlength):
        if (np_matrix[i][0], np_matrix[i][1]) in G.edges():
            for j in range(i):
                if np_matrix[j][0] == np_matrix[i][0] and \
                    np_matrix[j][1] == np_matrix[i][1]:
                    edgesrepeats[j] += 1
                    break   #add repeat number to the first edge
            G[np_matrix[i][0]][np_matrix[i][1]]['weight'] += float(np_matrix[i][2])
        else:
            G.add_edge(np_matrix[i][0], np_matrix[i][1], weight = float(np_matrix[i][2]))
    for i in range(matrixlength):
        if edgesrepeats[i] > 1:
            oldweight = G[np_matrix[i][0]][np_matrix[i][1]]['weight']
            G[np_matrix[i][0]][np_matrix[i][1]]['weight'] = oldweight / edgesrepeats[i]
    print("# of rows and columns in ml_matrix: %d %d" % (len(nodelist), len(nodelist)))
    print(nx.info(G))

    out_file = 'results/' + path.replace('.txt','_edges') + '.txt'
    print("writing %s to: %s" % (mlfile, out_file))
    with open(out_file, 'w') as out:
        if mlfile == "coactivation": #coactivation is undirected
            out.write('#node1\tnode2\tedge_weight\n')
        else:
            out.write('#tail\thead\tedge_weight\n')
        out.write('\n'.join(['%s\t%s\t%s' % (u,v,str(G.edge[u][v]['weight'])) for u,v in sorted(G.edges())]))
    return G

#BC-1 clustering coefficient for directed graph
def clustering_coefficient_directed_graph(G):
    """for directed graph, neibors include successors and predecessors"""
    cluster = {}
    numnodes = len(G.nodes())
    for node in G.nodes():
        succ = G.successors(node)
        pred = G.predecessors(node)
        neig = set(succ + pred)
        if len(neig) == 0 or len(neig) == 1:
            cluster[node] = 0
        else: 
            Gsub = nx.subgraph(G, list(neig))
            realedges = len(Gsub.edges())
            totaledges = len(neig) * (len(neig) - 1)
            cluster[node] = realedges * 1.0 / totaledges
    coeff = sum(cluster.values()) / len(cluster)
    return coeff

#BC-1 plot l and c for all graphs
def BC_plot_lc(lglist, cglist):
    """plot l and c for all six graphs"""
    yl, yc, xname = [], [], []
    x = range(1, len(lglist)+1)
    for i in lglist.keys():
        xname.append("  "+YEARS[i]+"   ")
        yl.append(lglist[i])
        yc.append(cglist[i])
    plt.title("Plot l and c for all six networks")
    plt.xlim([0, len(lglist)+1])
    plt.plot(x, yl, linestyle = 'None', marker = 's', label = 'l(G)')
    plt.plot(x, yc, linestyle = 'None', marker = 'o', label = 'c(G)')
    plt.legend()
    plt.xlabel(xname)
    plt.savefig('results/BC1_l_c.png') # save as png
    plt.clf()
    print()
    print("BC1_plot l and c for all six networks saved.")
    return

#BC-1 draw histogram for each network
def BC_plot_hist(G, mlfile):
    """plot histogram of l and c for all nodes in one network"""
    #calculate l and c value for each node
    lnodes = {}
    avelength = nx.shortest_path_length(G)
    for node, value in avelength.items():
        lnodes[node] = sum(value.values()) / len(value)

    if mlfile == 'coactivation':    #undirected graph
        cnodes = nx.clustering(G)
    else:
        cnodes = {}
        numnodes = len(G.nodes())
        for node in G.nodes():
            succ = G.successors(node)
            if len(succ) == 0 or len(succ) == 1:
                cnodes[node] = 0
            else: 
                Gsub = nx.subgraph(G, succ)
                realedges = len(Gsub.edges())
                totaledges = len(succ) * (len(succ) - 1)
                cnodes[node] = realedges * 1.0 / totaledges
    #plot histogram of l
    plt.hist(list(lnodes.values()))
    plt.title("l({}) for {} nodes".format(mlfile, len(lnodes)))
    plt.xlabel('l for each node')
    plt.ylabel('frequency')
    plt.savefig('results/BC1_hist_{}_l.png'.format(YEARS[mlfile]))
    plt.clf()
    print("Histogram for %s_l.png created." %YEARS[mlfile])

    #plot histogram of c
    plt.hist(list(cnodes.values()))
    plt.title("c({}) for {} nodes".format(mlfile, len(cnodes)))
    plt.xlabel('c for each node')
    plt.ylabel('frequency')
    plt.savefig('results/BC1_hist_{}_c.png'.format(YEARS[mlfile]))
    plt.clf()
    print("Histogram for %s_c.png created." %YEARS[mlfile])
    return

#BC-2 generate figure 1A of Cortical High-Density paper
def BC_figure1A_plot(Glist, lglist, filelist):
    """plot average path length against the denstiy of the network"""
    densitylist = {}
    for mlfile in filelist:
        G = Glist[mlfile]
        numedges = len(G.edges())
        numnodes = len(G.nodes())
        if mlfile == 'coactivation': #undirected graph
            densitylist[mlfile] = numedges * 2 / numnodes / (numnodes - 1)
        else:                      #derected graph
            densitylist[mlfile] = numedges / numnodes / (numnodes - 1)
        plt.plot(densitylist[mlfile], lglist[mlfile], linestyle = 'None', marker = 'o', label = mlfile+' '+YEARS[mlfile])
    
    #plt.ylim([1.0, 3.0])
    plt.title('Plot l and density')
    plt.xlabel('Graph density')
    plt.ylabel('Average path length')
    plt.legend()
    plt.savefig('results/BC2_density.png')
    #plt.clf()
    print("Figure 1A of Cortical High-Density paper created.")
    return densitylist

#BC-3 
def BC_figure1A_gray(G, densityG):
    """plot the gray area by randomly deleting edges"""
    densitytest = [x/10.0 for x in range(1, 7)]
    numden = len(densitytest)
    lave, lstd = [] * numden, [] * numden     #store the result for each density
    #G = parse_high_density_file(mlfile, path)
    alledges = G.edges()
    numedges = G.size()
    numnodes = G.order()
    maxedges = numnodes * (numnodes - 1)
    lave = []
    lstd = []

    for density in densitytest:
        removeedges = math.floor(maxedges * (densityG - density))
        ld = []
        for n in range(100):   #repeat 100 times
            Gremove = nx.DiGraph()
            Gremove.add_edges_from(alledges)
            removeloc = random.sample(range(numedges), removeedges)
            for i in removeloc:
                u = alledges[i][0]
                v = alledges[i][1]
                Gremove.remove_edge(u, v)

            #default average_shortest_path_lenth function cannot handle unconnected graph
            lnodes = {}
            avelength = nx.shortest_path_length(Gremove)
            for node, value in avelength.items():
                lnodes[node] = sum(value.values()) / len(value)
            l = sum(lnodes.values()) / len(lnodes)
            ld.append(l)
        npld = np.array(ld)
        lave.append(npld.mean())
        lstd.append(npld.std())

    #draw graph
    plt.errorbar(densitytest, lave, yerr = lstd, marker = 's', markersize = 2, label = 'gray area')
    plt.title('Plot l and density with gray area')
    plt.legend()
    plt.savefig('results/BC3_gray_area.png')
    plt.clf()
    print("Figure 1A with gray area created.")
    return lave, lstd