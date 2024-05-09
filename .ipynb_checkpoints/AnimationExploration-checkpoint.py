import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random

# Create a random graph
G = nx.erdos_renyi_graph(10, 0.3)

# Initialize the figure and set up the plot
fig, ax = plt.subplots()
pos = nx.spring_layout(G)  # Define the layout

# Function to update the node colors in each frame
def update(frame):
    colors = []  # Initialize empty list for node colors
    for node in G.nodes():
        # Generate random RGB values for node color
        random_color = (random.random(), random.random(), random.random())
        colors.append(random_color)
    
    # Clear the previous plot and update node colors
    ax.clear()
    nx.draw(G, pos, with_labels=True, node_color=colors, ax=ax)

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=10, interval=1000, repeat=True)

plt.show()