
import networkx as nx


def get_node_data(f_graph):
    g_name_to_index = {}
    g_index_to_name = {}
    nodes_count = 0
    index = 0
    for i, row in enumerate(f_graph):
        if i < 4:
            continue
        if i == 4:
            nodes_count = int(row) + i
            continue
        if i == nodes_count + 1:
            print row
            break
        g_name_to_index[row] = index
        g_index_to_name[index] = row
        index += 1
    return g_name_to_index, g_index_to_name

def clean_node(node):
    return node[2:-2]


if __name__ == "__main__":
    graph_name = 'syeast0_del_nodes_100_modify_0.05_iter_0'
    # f_graph1 = open('output_graphs/{}/{}.gw'.format(graph_name, graph_name))
    f_graph1 = open('output_alignment/graphs/{}.gw'.format(graph_name))
    g1_name_to_index, g1_index_to_name = get_node_data(f_graph1)

    # f_graph2 = open('networks/yeast/yeast.gw')
    # g2_name_to_index, g2_index_to_name = get_node_data(f_graph2)

    ouput_file = open('%s.align' % graph_name, 'wb')
    for i in range(len(g1_index_to_name)):
        ouput_file.write(clean_node(g1_index_to_name[i].rstrip()) + '\t' + clean_node(g1_index_to_name[i].rstrip()) + '\n')
    # graph2 = nx.read_leda('networks/yeast/yeast.gw')
    pass
