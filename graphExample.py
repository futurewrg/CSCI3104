#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 15:29:04 2017

@author: rhoenigman
"""

import networkx as nx
import matplotlib.pyplot as plt  
#create a directed graph
G = nx.DiGraph()

#adding an edge also adds the node
G.add_edge('s', 'p1', weight=1.0)
G.add_edge('s', 'p2', weight=1.0)
G.add_edge('s', 'p3', weight=1.0)

G.add_edge('p1', 'i1', weight=1.0)
G.add_edge('p1', 'i2', weight=1.0)

G.add_edge('p2', 'i2', weight=1.0)
G.add_edge('p2', 'i3', weight=1.0)

G.add_edge('p3', 'i1', weight=1.0)
G.add_edge('p3', 'i3', weight=1.0)

G.add_edge('i1', 't', weight=1.0)
G.add_edge('i2', 't', weight=1.0)
G.add_edge('i3', 't', weight=1.0)

#each edge has a weight of 1. The shortest path is the fewest edges.
#Use this to verify that your graph built correctly.
#
t = nx.all_simple_paths(G, 'Spider', 'Fly')

nx.draw(G, with_labels=True)  
plt.show()




