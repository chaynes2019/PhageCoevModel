from pyvis.network import Network
import networkx as nx


nt = Network(notebook=False)
nt.add_node(1)
nt.add_node(2)
nt.show("example.html")