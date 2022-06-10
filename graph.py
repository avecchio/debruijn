import random
import toyplot
from Bio import SeqIO

sequence = "10101101111000101101010110011"
#with open("sequence.fasta") as handle:
#    for record in SeqIO.parse(handle, "fasta"):
#        sequence = str(record.seq)

def get_kmer_count_from_sequence(sequence, k=3, cyclic=True):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence.
    """
    # dict to store kmers
    kmers = {}
    
    # count how many times each occurred in this sequence (treated as cyclic)
    for i in range(0, len(sequence)):
        kmer = sequence[i:i + k]
        
        # for cyclic sequence get kmers that wrap from end to beginning
        length = len(kmer)
        if cyclic:
            if len(kmer) != k:
                kmer += sequence[:(k - length)]
        
        # if not cyclic then skip kmers at end of sequence
        else:
            if len(kmer) != k:
                continue
        
        # count occurrence of this kmer in sequence
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    
    return kmers

def get_debruijn_edges_from_kmers(kmers):
    """
    Every possible (k-1)mer (n-1 suffix and prefix of kmers) is assigned
    to a node, and we connect one node to another if the (k-1)mer overlaps 
    another. Nodes are (k-1)mers, edges are kmers.
    """
    # store edges as tuples in a set
    edges = set()
    
    # compare each (k-1)mer
    for k1 in kmers:
        for k2 in kmers:
            if k1 != k2:            
                # if they overlap then add to edges
                if k1[1:] == k2[:-1]:
                    edges.add((k1[:-1], k2[:-1]))
                if k1[:-1] == k2[1:]:
                    edges.add((k2[:-1], k1[:-1]))

    return edges

import toyplot.png
import networkx as nx
import matplotlib.pyplot as plt

def plot_debruijn_graph(edges, width=500, height=500):

    G = nx.DiGraph()
    G.add_edges_from(edges)
    #    [('A', 'B'), ('A', 'C'), ('D', 'B'), ('E', 'C'), ('E', 'F'),
    #    ('B', 'H'), ('B', 'G'), ('B', 'F'), ('C', 'G')])

    #val_map = {'A': 1.0,
    #        'D': 0.5714285714285714,
    #        'H': 0.0}

    #values = [val_map.get(node, 0.25) for node in G.nodes()]

    # Specify the edges you want here
    #red_edges = [('A', 'C'), ('E', 'C')]
    #edge_colours = ['black' if not edge in red_edges else 'red'
    #                for edge in G.edges()]
    #black_edges = [edge for edge in G.edges() if edge not in red_edges]

    # Need to create a layout when doing
    # separate calls to draw nodes and edges
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), 
                        #node_color = values, 
                        node_size = 500)
    nx.draw_networkx_labels(G, pos)
    #nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', arrows=True)
    nx.draw_networkx_edges(G, pos, edgelist=edges, arrows=False)
    plt.savefig("graph.png")
    '''"returns a toyplot graph from an input of edges"
    graph = toyplot.graph(
        [i[0] for i in edges],
        [i[1] for i in edges],
        width=width,
        height=height,
        tmarker=">", 
        vsize=25,
        vstyle={"stroke": "black", "stroke-width": 2, "fill": "none"},
        vlstyle={"font-size": "11px"},
        estyle={"stroke": "black", "stroke-width": 2},
        layout=toyplot.layout.FruchtermanReingold(edges=toyplot.layout.CurvedEdges()))
    toyplot.png.render(graph, "grpah.png")
    return graph
    '''

kmers = get_kmer_count_from_sequence(sequence, k=3)
print(kmers)
edges = get_debruijn_edges_from_kmers(kmers)
plot_debruijn_graph(list(edges))