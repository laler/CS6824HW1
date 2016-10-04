###################################################################
#CS6824 Homework1
#Zhen Guo
#10/2/16
#Draw Figure 1 of Collective dynamics of 'small-world' networks.
#Draw Figure 1A of Cortical High-Density Counterstream Architecture. 
###################################################################

import WSGraph
import BrainConnectomes
import networkx as nx
import numpy as np
import pylab as plt
import time

MLFILES = {
    "primate_cortex": "./1991-cerebral-cortex-felleman-primate-cerebral-cortex-fv30.mat",
    "macaque_cortex": "./1993-Proc-Royal-Society-organization-neural-systems-macaque71.mat",
    "cat_cortex": "./1995-journal-neuroscience-connectivity-cerebral-cortex-cat.mat",
    "macaque_functional": "./2007-pnas-honey-network-structure-functional-connectivity-macaque47.mat",
    "coactivation": "./2013-PNAS-Crossley-Cognitive-relevance-coactivation-matrix.mat",
    "high_density": "./2012-cerebral-cortex-markov-weighted-directed-interareal-macaque.txt"
}

#WS-3 run programs using two sets of p
nlist = [100, 1000, 5000]
k = 20
p1 = [0.05 * k for k in range(1, 21)]
p2 = WSGraph.p_range(k)
f = open('results/WS3_p1.out', 'w')
WSGraph.test_p_n(nlist, k, p1, f)
f.close()
f = open('results/WS3_p2.out', 'w')
WSGraph.test_p_n(nlist, k, p2, f)
f.close()

#WS-4 generate Fig.1 in small-world paper
for n in nlist:
    WSGraph.SW_figure1_plot(n, k, p2)

#BC-1 parse data from six network files
Glist = {}  #store six graphs
for mlfile in MLFILES:
    if mlfile != "high_density":
        Glist[mlfile] = BrainConnectomes.parse_mat_files(mlfile, MLFILES[mlfile])
    else: #parse the last graph
        Glist[mlfile] = BrainConnectomes.parse_high_density_file(mlfile, MLFILES[mlfile])

#BC-1 plot l and c for all graphs
lglist, cglist = {}, {}
for mlfile in MLFILES:
    lglist[mlfile] = WSGraph.average_shortest_path_length(Glist[mlfile])
    if mlfile == 'coactivation':  #undirected graph
        cglist[mlfile] = WSGraph.clustering_coefficient(Glist[mlfile])
    else:                         #directed graph
        cglist[mlfile] = BrainConnectomes.clustering_coefficient_directed_graph(Glist[mlfile])
BrainConnectomes.BC_plot_lc(lglist, cglist)
               
#BC-1 plot histogram for each graph
for mlfile in MLFILES:
    BrainConnectomes.BC_plot_hist(Glist[mlfile], mlfile)

#BC-2 generate Fig.1A in Cortical high-density paper
densitylist = BrainConnectomes.BC_figure1A_plot(Glist, lglist, list(MLFILES.keys()))
f = open('results/BC2_l_c_d.out','w')
f.write("average shortest path length: \n")
f.write(str(lglist))
f.write("\naverage clustering coefficient: \n")
f.write(str(cglist))
f.write("\ngraph density: \n")
f.write(str(densitylist))
f.close()
print("Writing l, c, d values for six graphs completed.")

#BC-3 plot gray area for Fig.1A
mlfile = "high_density"
lave, lstd = BrainConnectomes.BC_figure1A_gray(Glist[mlfile], densitylist[mlfile])

