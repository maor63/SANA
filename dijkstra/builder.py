from graph import Graph
import pickle
import numpy as np
import lzma
"""
Don't try to replace the pickle function with the numpy save/load.
It's 1) a waste of time and 2) less space efficient
"""
"""
build_graph() takes a file in the following format:
Node 1 Node 2
and returns a graph with every edge Node 1 <-> Node 2
"""


def convert_graph_int(file):
	nodes = {}
	count = 0
	int_graph = open("int_graph.txt", "w")
	for line in open(file):
		n1, n2 = line.split()
		if n1 not in nodes:
			nodes[n1] = count
			count += 1
		if n2 not in nodes:
			nodes[n2] = count
			count += 1
		int_graph.write(f"{nodes[n1]} {nodes[n2]}\n")

	dict_file = open("dict.txt", 'w')
	for node, value in nodes.items():
		dict_file.write(f"{value} {node}\n")

	int_graph.close()
	dict_file.close()
		

def build_graph(file):
	count = 0
	g = Graph()
	for line in open(file):
		n1, n2 = line.split()
		if n1 not in g.indexes:
			g.indexes[n1] = count
			g.nodes[count] = n1
			count += 1
		if n2 not in g.indexes:
			g.indexes[n2] = count
			g.nodes[count] = n2
			count += 1
		g.add_edge(g.indexes[n1], g.indexes[n2])
	return g

"""
get_sim() provides a wrapper to speed up the construction of the sims matrix.
It stores the matrix in a pickle to minimize calls to build_sim. Future calls
to this function will just load the pickle, which is faster for large datasets
than calling build_sims().
"""
def get_sim(file, graph1, graph2, pickle_name = ".sim.pickle"):
    """
    build_sim takes a file and builds a [yeast_size X human_size]
    array of similarity values. This function takes a long time to run
    for large datasets so it should be avoided when possible.
    runtime: O(|Y| * |H|)
    """
    def build_sim(file, graph1, graph2):
        similarity = np.zeros((len(graph1), len(graph2)))
        with lzma.open(file, mode = 'rt') as f:
            for line in f:
                node_y, node_h, sim_val = line.strip().split() # 8 seconds
                try:
                    node_y, node_h, sim_val = graph1.indexes[node_y], graph2.indexes[node_h], float(sim_val) #17 seconds
                    similarity[node_y][node_h] = sim_val
                except:
                    pass
        return similarity
    
    if ".sim.pickle" == pickle_name:
        pickle_name = graph1.name + graph2.name + pickle_name
    try:
        with open(pickle_name,'rb') as f:
            return pickle.load(f)
    except FileNotFoundError as e:
        sims = build_sim(file, graph1, graph2)
        with open(pickle_name,'wb') as f:
            pickle.dump(sims,f)
        return sims
"""
functions below are not working, stubs for future releases
"""
def read_dict(file):
    d = dict()
    with open(file, 'r') as f:
        for line in f:
            row = line.strip().split(" ")
            d[int(row[0])] = row[1]
    return d


