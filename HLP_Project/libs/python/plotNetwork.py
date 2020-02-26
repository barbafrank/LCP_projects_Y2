import networkx as nx
import scipy
import numpy as np
import random
from matplotlib import pyplot as plt
import matplotlib as mpl

def plotNetworkClusters(A, clusters_vector, node_size=10, figsize=(15, 7),
                                        draw_edges=True, pos=None, colors=None, ax=None, cbar=True, return_vals=False):

    clusters = {k: v for k, v in enumerate(clusters_vector)}
    comms = set(clusters.values())
    size = len(comms)

    #Generating a random set of colors (one for each cluster)
    if colors is None:
        colors = generateRandomColors(size)

    # test creating the cmap for the colorbar
    cmap = mpl.colors.ListedColormap(colors)
    norm = mpl.colors.BoundaryNorm(np.arange(-1, size)-0.5, cmap.N)
    #print(cmap)

    #Creating graph
    G = nx.from_numpy_matrix(A)

    #Plotting
    #plt.figure(num=None, figsize=figsize, facecolor='w', edgecolor='k')
    if ax is None:
        fig, ax = plt.subplots(1, 1, num=None, figsize=figsize, facecolor='w', edgecolor='k')
    if cbar:
        bar_ax, _ = mpl.colorbar.make_axes(ax, 'bottom')

    if pos is None:
        pos = nx.spring_layout(G)

    for i, com in enumerate(comms):
        list_nodes = [nodes for nodes in clusters.keys() if clusters[nodes] == com]
        nx.draw_networkx_nodes(G, pos, list_nodes, node_size, node_color = colors[i], ax=ax)
        #nx.draw_networkx_nodes(G, pos, list_nodes, node_size, node_color = rgb2hex(colors[i][:3]))

    if(draw_edges):
        nx.draw_networkx_edges(G, pos, alpha=0.5, ax=ax)

    # set the colorbar
    if cbar:
        cb2 = mpl.colorbar.ColorbarBase(bar_ax, cmap=cmap,
                                        norm=norm,
                                        ticks=np.arange(-1, size),
                                        spacing='proportional',
                                        orientation='horizontal')
        bar_ax.set_xticklabels(list(comms))

    if return_vals:
        return pos, colors

def generateRandomColors(number_of_colors):
    return ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]
