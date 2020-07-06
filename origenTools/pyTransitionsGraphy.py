import numpy as np
from pyOrigen import ORIGENData
import pandas as pd
from periodictable import elements
import sys
import matplotlib.pyplot as plt
import graph_tool.all as gt

diagDecayFile = sys.argv[1]
diagRxFile = sys.argv[2]
offDiagRxFile = sys.argv[3]

data = ORIGENData(diagDecayFile, diagRxFile, offDiagRxFile)
ele = elements.symbol('Xe').number

XeIDs = data.getNuclideOrigenIDs('Xe', addG1=False)
IodineIDs = data.getNuclideOrigenIDs('I', addG1=False)
"""
LiIDs = data.getNuclideOrigenIDs('Li', 7, addG1=False)
FIDs = data.getNuclideOrigenIDs('F', 19, addG1=False)
BeIDs = data.getNuclideOrigenIDs('Be', 9, addG1=False)
ZrIDs = data.getNuclideOrigenIDs('Zr', 90, addG1=False)
KrIDs = data.getNuclideOrigenIDs('Kr', addG1=False)
NbIDs = data.getNuclideOrigenIDs('Nb', addG1=False)
MoIDs = data.getNuclideOrigenIDs('Mo', addG1=False)
CdIDs = data.getNuclideOrigenIDs('Cd', addG1=False)
NdIDs = data.getNuclideOrigenIDs('Nd', addG1=False)
SmIDs = data.getNuclideOrigenIDs('Sm', addG1=False)
"""

"""
OrigenIDs = [XeIDs, IodineIDs, LiIDs, FIDs, BeIDs, ZrIDs, 
    KrIDs, NbIDs, MoIDs, CdIDs, NdIDs, SmIDs]

"""
OrigenIDs = [XeIDs, IodineIDs]
#for atomicNumber in range(40,57):
#   IDs = data.getNuclideOrigenIDs(atomicNumber, addG1=False)
#   OrigenIDs.append(IDs)

G = gt.Graph(directed=True)
vOrigenID = G.new_vertex_property("int")
vAtomicNumber = G.new_vertex_property("int")
#eReactionRate = G.new_edge_property("double")
name_v_map = {}

for IDs in OrigenIDs:
    for ID in IDs:
        if ID not in name_v_map:
            myv = G.add_vertex()
            vOrigenID[myv] = ID
            vAtomicNumber[myv] = data.getElement(ID)
            name_v_map[ID] = myv
        parents = data.getReactionParents(ID)
        parents = data.keepActinide(parents, 92, 235)
        for parent in parents:
            if parent not in name_v_map:
                pv = G.add_vertex()
                vOrigenID[pv] = parent
                vAtomicNumber[pv] = data.getElement(parent)
                name_v_map[parent] = pv
        daughters = data.getReactionDaughters(ID)
        daughters = data.keepActinide(daughters, 92, 235)
        for daughter in daughters:
            if daughter not in name_v_map:
                dv = G.add_vertex()
                vOrigenID[dv] = daughter
                vAtomicNumber[dv] = data.getElement(daughter)
                name_v_map[daughter] = dv

for IDs in OrigenIDs:
    for ID in IDs:
        myv = name_v_map[ID]
        parents = data.getReactionParents(ID)
        parents = data.keepActinide(parents, 92, 235)
        for parent in parents:
            pv = name_v_map[parent]
            G.add_edge(pv, myv)
        daughters = data.getReactionDaughters(ID)
        daughters = data.keepActinide(daughters, 92, 235)
        for daughter in daughters:
            dv = name_v_map[daughter]
            G.add_edge(myv, dv)
SIZE = 1500
V_SIZE = SIZE / 70
E_PWIDTH = V_SIZE/12
#deg = G.degree_property_map("in")
#deg.a = V_SIZE + 5.*deg.a
G.purge_vertices()
state = gt.minimize_nested_blockmodel_dl(G, deg_corr=True)
t = gt.get_hierarchy_tree(state)[0]
tpos = pos = gt.radial_tree_layout(t, t.vertex(t.num_vertices() - 1), weighted=True)
cts = gt.get_hierarchy_control_points(G, t, tpos)
pos = G.own_property(tpos)
b = state.levels[0].b
shape = b.copy()
shape.a %= 14


gt.graph_draw(G,
    pos = pos,
    output_size = (SIZE,SIZE),
    vertex_size = V_SIZE,
    edge_control_points=cts,
    edge_pen_width = E_PWIDTH,
    #vertex_fill_color = vAtomicNumber,
    #vertex_text = vAtomicNumber
    #output="dopeFPpic.pdf"
)
