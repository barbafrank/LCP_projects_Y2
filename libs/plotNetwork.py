import networkx as nx
import scipy
import numpy as np
import random
from matplotlib import pyplot as plt

def plotNetworkClusters(A, clusters_vector, node_size=10, figsize=(10, 5), draw_edges=True):
    
    size = len(clusters_vector)
    clusters = {k: v for k, v in enumerate(clusters_vector)}
    
    #Generating a random set of colors (one for each cluster)
    colors = generateRandomColors(size+1)
    
    #Creating graph
    G = nx.from_numpy_matrix(A)
    
    #Plotting
    plt.figure(num=None, figsize=figsize, facecolor='w', edgecolor='k')
    pos = nx.spring_layout(G)
    
    for com in set(clusters.values()):
        list_nodes = [nodes for nodes in clusters.keys() if clusters[nodes] == com]
        nx.draw_networkx_nodes(G, pos, list_nodes, node_size, node_color = colors[com+1])
    
    if(draw_edges):
        nx.draw_networkx_edges(G, pos, alpha=0.5)
    
    plt.show()
    
def generateRandomColors(number_of_colors):
    return ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]