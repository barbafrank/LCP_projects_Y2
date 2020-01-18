import matplotlib.pyplot as plt
from libs.approximateSimrank import approximateSimrank
from libs.edgelistParser import edgelistParser


# load the adjacency matrix in the expected format
A, N, D = edgelistParser("datasets/facebook/0.edges")

# set source and teleportation rate
c = .85
v = 10
epsilon = 1e-5

# compute the approximateSimrank
p = approximateSimrank(A, N, D, v, c, epsilon)

fig, axs = plt.subplots(1, 1, figsize=(7,7))
axs.plot(p, '.')
plt.show()
