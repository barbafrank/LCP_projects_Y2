import networkx as nx
import scipy
import numpy as np
import random
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import rgb2hex

def plotNetworkClusters(A, clusters_vector, node_size=10, figsize=(10, 5), draw_edges=True, pos=None):

    size = len(clusters_vector)
    clusters = {k: v for k, v in enumerate(clusters_vector)}

    #Generating a random set of colors (one for each cluster)
    #colors = generateRandomColors(size+1)
    x = np.arange(size+1)
    colors = cm.get_cmap('tab20')(x)

    #Creating graph
    G = nx.from_numpy_matrix(A)

    #Plotting
    plt.figure(num=None, figsize=figsize, facecolor='w', edgecolor='k')
    if pos is None:
        pos = nx.spring_layout(G)

    for i, com in enumerate(set(clusters.values())):
        list_nodes = [nodes for nodes in clusters.keys() if clusters[nodes] == com]
        nx.draw_networkx_nodes(G, pos, list_nodes, node_size, node_color = rgb2hex(colors[i][:3]))

    if(draw_edges):
        nx.draw_networkx_edges(G, pos, alpha=0.5)

    plt.show()
    return pos

def generateRandomColors(number_of_colors):
    return ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]
