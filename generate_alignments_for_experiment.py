from __future__ import print_function

from collections import Counter

import matplotlib.pyplot as plt
import timeit
from multiprocessing import Pool, Process
import networkx as nx
import random
import pandas as pd
import os
import numpy as np

from eval_true_alignment import get_node_data


def write_leda_to_file(G, graph_name, path):
    header = """LEDA.GRAPH\nstring\nlong\n-2\n"""
    nodes = G.nodes()
    node_index = {n: i + 1 for i, n in enumerate(nodes)}
    edges = G.edges()
    if not os.path.exists(path):
        os.makedirs(path)

    output_f = open(os.path.join(path, '{0}.gw'.format(graph_name)), 'wb')
    output_f.write(header)
    rows = []
    rows.append('{0}'.format(len(nodes)))
    for node in nodes:
        rows.append('|{' + str(node) + '}|')
    rows.append('{0}'.format(len(edges)))
    for v, u in edges:
        m = node_index[v]
        j = node_index[u]
        rows.append('{0} {1} 0 '.format(m, j) + '|{0}|')

    output_f.write('\n'.join(rows))
    output_f.write('\n')
    output_f.close()


def get_nc(alignment, true_alignment):
    correct = sum(1 for v in alignment.index if alignment[v] == true_alignment[v])
    return correct * 1.0 / len(alignment)


def update_Ea_change(v, u_old, u_new, G1, G2, A):
    res = 0
    for n in G1.neighbors(v):
        res -= 1 if G2.has_edge(u_old, A[n]) else 0
        res += 1 if G2.has_edge(u_new, A[n]) else 0
    return res


def update_Ea_hat(u_old, u_new, G2, set_free_nodes):
    res = 0
    for n in G2.neighbors(u_old):
        res -= 1 if n not in set_free_nodes else 0
    for n in G2.neighbors(u_new):
        res += 1 if n not in set_free_nodes else 0

    res -= 1 if G2.has_edge(u_old, u_new) else 0
    return res


def add_random_edges(G, k):
    assert isinstance(G, nx.Graph)
    edges_to_add = random.sample(list(nx.non_edges(G)), k)
    G.add_edges_from(edges_to_add)
    return G


def remove_random_edges(G, k):
    assert isinstance(G, nx.Graph)
    edges_to_remove = random.sample(list(G.edges), k)
    G.remove_edges_from(edges_to_remove)
    return G


def get_wattz_graphs():
    n = 1000
    WS = nx.watts_strogatz_graph(n, 18, 0.5)
    WS.graph['name'] = 'ws_1000_18_0.5'
    WS_bal = nx.watts_strogatz_graph(n, 500, 0.5)
    WS_bal.graph['name'] = 'ws_1000_500_0.5'
    WS_dense = nx.watts_strogatz_graph(n, 980, 0.5)
    WS_dense.graph['name'] = 'ws_1000_982_0.5'
    G2s = [WS, WS_bal, WS_dense]
    return G2s


def get_erdos_graphs():
    n = 1000
    ER = nx.gnm_random_graph(n, 8500)
    ER.graph['name'] = 'er_1000_8500'
    ER_bal = nx.gnm_random_graph(n, 250000)
    ER_bal.graph['name'] = 'er_1000_250000'
    ER_dense = nx.gnm_random_graph(n, 491000)
    ER_dense.graph['name'] = 'er_1000_491000'
    G2s = [ER, ER_bal, ER_dense]
    return G2s


def get_barabasi_graphs():
    n = 1000
    BA = nx.barabasi_albert_graph(n, 9)
    BA.graph['name'] = 'ba_1000_9'
    BA_bal = nx.barabasi_albert_graph(n, 500)
    BA_bal.graph['name'] = 'ba_1000_500'
    BA_compliment = nx.complement(nx.barabasi_albert_graph(n, 9))
    BA_compliment.graph['name'] = 'ba_1000_9_compliment'
    G2s = [BA, BA_bal, BA_compliment]
    return G2s


def syeast0_graphs():
    G2_name = 'syeast0'
    Y0 = nx.read_leda('output_graphs/{0}/{0}.gw'.format(G2_name))
    Y0.graph['name'] = G2_name
    Y0_complement = nx.complement(Y0)
    Y0_complement.graph['name'] = G2_name + "_complement"
    write_leda_to_file(Y0_complement, Y0_complement.graph['name'], 'output_alignment/graphs')
    Y0_V = len(Y0.nodes)
    edges_count_to_add = int(Y0_V * (Y0_V - 1) / 4.0 - len(Y0.edges))
    Y0_balanced = add_random_edges(nx.Graph(Y0), edges_count_to_add)
    Y0_balanced.graph['name'] = 'syeast0_balanced'
    write_leda_to_file(Y0_balanced, Y0_balanced.graph['name'], 'output_alignment/graphs')
    G2s = [Y0, Y0_complement, Y0_balanced]
    return G2s


def generate_alignments(G2s, edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, shuffles,
                        output_file_name='alignmnets_y0_com_bal'):
    iter_counter = 1
    for i in range(full_iterations):
        for j, G2 in enumerate(G2s):
            for nodes_to_delete in nodes_to_delete_list:
                for edge_frac_To_modify in edge_frac_To_modify_list:
                    start = timeit.default_timer()
                    G2_name = G2.graph['name']
                    print('full iteration {}/{}, using G2 {}, nodes_to_delete {}, edge_to_delete {} {}/{}'.format(
                        str(i + 1), full_iterations, G2_name, nodes_to_delete, edge_frac_To_modify, iter_counter,
                        str(len(G2s) * len(nodes_to_delete_list) * len(edge_frac_To_modify_list) * full_iterations)))
                    iter_counter += 1
                    G1 = generate_G1(G2, nodes_to_delete, edge_frac_To_modify)
                    G1_name = '{}_del_nodes_{}_modify_{}_iter_{}'.format(G2.graph['name'], nodes_to_delete,
                                                                         edge_frac_To_modify, i)
                    write_leda_to_file(G1, G1_name, 'output_alignment/graphs')

                    # generate(G1, G1_name, G2, G2_name, edge_frac_To_modify, i, nodes_to_delete, output_file_name,
                    #          shuffles, start)


def generate(G1, G1_name, G2, G2_name, edge_frac_To_modify, i, nodes_to_delete, output_file_name, shuffles, start):
    V1 = list(G1.nodes)
    V2 = list(G2.nodes)
    E2 = len(list(G2.edges))
    true_A = pd.Series(data=V1, index=V1)
    E1, Ea, Ea_hat = get_align_edges(G1, G2, true_A)
    Ea_hat_star = Ea_hat
    Ea_star = Ea
    omega = (len(V1) * (len(V1) - 1)) / 2.0
    try:
        nc = get_nc(A, true_A)
        measures = get_measures(E1, Ea, Ea_hat, omega)
        acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measures

        row = (
            G1_name, len(V1), E1, G2_name, len(V2), E2, omega, Ea_hat_star, Ea, Ea_hat, '', 0, '', 0,
            '', 0,
            nc, tp, fp, tn, fn, ec, ics, s3, f1, acc, bm, mk, mcc, f01, f033, f3, f10)
        rows = [row]
    except Exception as e:
        rows = []
    original_Ea = Ea
    original_Ea_hat = Ea_hat
    set_free_nodes = set(G2.nodes) - set(G1.nodes)
    for shuffle in range(shuffles):
        print('\r shuffle iteration {}/{}'.format(shuffle, shuffles), end='')
        A = pd.Series(true_A.copy())
        random.shuffle(V1)
        G1_name = '{}_del_nodes_{}_modify_{}_iter_{}_shuffle_{}'.format(G2.graph['name'],
                                                                        nodes_to_delete,
                                                                        edge_frac_To_modify, i,
                                                                        shuffle)
        Ea = original_Ea
        Ea_hat = original_Ea_hat

        for v in V1:
            u = A[v]
            u2 = random.choice(list(set_free_nodes))

            ea_add = update_Ea_change(v, u, u2, G1, G2, A)
            ea_hat_add = update_Ea_hat(u, u2, G2, set_free_nodes)

            Ea += ea_add
            Ea_hat += ea_hat_add

            change(A, v, u2, set_free_nodes)
            try:
                nc = get_nc(A, true_A)
                measuers = get_measures(E1, Ea, Ea_hat, omega)
                acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measuers

                row = (
                    G1_name, len(V1), E1, G2_name, len(V2), E2, omega, Ea_hat_star, Ea, Ea_hat, v,
                    G1.degree[v], u,
                    G2.degree[u], u2, G2.degree[u2], nc, tp, fp, tn, fn, ec, ics, s3, f1, acc, bm, mk,
                    mcc,
                    f01,
                    f033, f3, f10)
                rows.append(row)
            except Exception as e:
                print(e)
    stop = timeit.default_timer()
    print('\r shuffle iteration {}/{} time: {}sec'.format(shuffles, shuffles, stop - start), end='')
    print()
    df = pd.DataFrame(rows,
                      columns=['G1', 'V1', 'E1', 'G2', 'V2', 'E2', 'omega', 'Ea_hat*', 'Ea', 'Ea_hat',
                               'G1 changed node', 'Degree of changed node', 'G2 old node',
                               'Degree of old node', 'G2 new node', 'Degree of new node', 'NC', 'TP',
                               'FP', 'TN', 'FN', 'EC', 'ICS', 'S3', 'F1', 'ACC', 'BM', 'MK', 'MCC',
                               'F0.1', 'F0.33', 'F3', 'F10'])
    output_path = 'output_alignment/%s.csv' % output_file_name
    if os.path.isfile(output_path):
        df.to_csv(output_path, mode='a', header=None)
    else:
        df.to_csv(output_path)


def change(A, v, u2, set_free_nodes):
    u = A[v]
    A[v] = u2
    set_free_nodes.remove(u2)
    set_free_nodes.add(u)


def get_measures(E1, Ea, Ea_hat, omega):
    ec, f1, ics, s3, bm, mk, mcc = compute_measures(E1, Ea, Ea_hat, omega)
    fn, fp, tn, tp = get_confution_matrix(E1, Ea, Ea_hat, omega)
    acc = (tp + tn) / omega * 1.0
    f01 = get_f_beta(E1, Ea, Ea_hat, 0.1)
    f033 = get_f_beta(E1, Ea, Ea_hat, 1.0 / 3)
    f3 = get_f_beta(E1, Ea, Ea_hat, 3)
    f10 = get_f_beta(E1, Ea, Ea_hat, 10)
    measures = [acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp]
    measure_names = ['acc', 'bm', 'ec', 'f01', 'f033', 'f1', 'f10', 'f3', 'fn', 'fp', 'ics', 'mcc', 'mk', 's3',
                     'tn', 'tp']
    return pd.Series(index=measure_names, data=measures)


def get_f_beta(E1, Ea, Ea_hat, beta):
    f_beta = ((1 + beta ** 2) * Ea) / float(E1 + beta ** 2 * Ea_hat)
    return f_beta


def get_confution_matrix(E1, Ea, Ea_hat, omega):
    tp = Ea
    fp = E1 - Ea
    fn = Ea_hat - Ea
    tn = omega - (tp + fp + fn)
    return fn, fp, tn, tp


def compute_measures(E1, Ea, Ea_hat, omega):
    ec = Ea / E1
    ics = Ea / Ea_hat
    s3 = Ea / (E1 + Ea_hat - Ea)
    f1 = (2 * Ea) / (E1 + Ea_hat)
    bm = (omega * Ea - E1 * Ea_hat + 1.0) / ((omega - Ea_hat) * Ea_hat + 1.0)
    mk = (omega * Ea - E1 * Ea_hat + 1.0) / ((omega - E1) * E1 + 1.0)
    mcc = (omega * Ea - E1 * Ea_hat + 1.0) / ((E1 * Ea_hat * (omega - E1) * (omega - Ea_hat)) ** 0.5 + 1.0)
    # mcc_geometric = (bm + mk) / 2.0
    return ec, f1, ics, s3, bm, mk, mcc


def generate_G1(Y0, nodes_to_delete, edges_to_modify):
    G1 = nx.Graph(Y0)
    v1 = list(G1.nodes)
    G1.remove_nodes_from(random.sample(v1, nodes_to_delete))

    total_edges_to_change = int(len(v1) * edges_to_modify)
    edges_to_add = random.randint(0, total_edges_to_change)
    edges_to_remove = total_edges_to_change - edges_to_add
    add_random_edges(G1, edges_to_add)
    remove_random_edges(G1, edges_to_remove)
    return G1


def get_align_edges(G1, G2, A):
    Ea_hat = len(G2.subgraph(A.values).edges) * 1.0
    Ea = sum([1 if G2.has_edge(A[v1], A[v2]) else 0 for v1, v2 in G1.edges]) * 1.0
    E1 = len(G1.edges) * 1.0
    return E1, Ea, Ea_hat


def update_Ea_swap(v1, v2, u1, u2, G1, G2, A):
    res = 0
    for n in G1.neighbors(v1):
        res -= 1 if G2.has_edge(u1, A[n]) else 0
        res += 1 if G2.has_edge(u2, A[n]) else 0

    for n in G1.neighbors(v2):
        res -= 1 if G2.has_edge(u2, A[n]) else 0
        res += 1 if G2.has_edge(u1, A[n]) else 0

    res += 2 if G2.has_edge(u1, u2) and G1.has_edge(v1, v2) else 0
    return res


def generate_alignments_for_objective(G2s, edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, objective,
                                      output_file_name='alignmnets_y0_com_bal'):
    iter_counter = 1
    for i in range(full_iterations):
        for j, G2 in enumerate(G2s):
            for nodes_to_delete in nodes_to_delete_list:
                for edge_frac_To_modify in edge_frac_To_modify_list:

                    G2_name = G2.graph['name']
                    print('full iteration {}/{}, using G2 {}, nodes_to_delete {}, edge_to_delete {} {}/{}'.format(
                        str(i + 1), full_iterations, G2_name, nodes_to_delete, edge_frac_To_modify, iter_counter,
                        str(len(G2s) * len(nodes_to_delete_list) * len(edge_frac_To_modify_list) * full_iterations)))
                    iter_counter += 1
                    G1 = generate_G1(G2, nodes_to_delete, edge_frac_To_modify)
                    G1_name = '{}_del_nodes_{}_modify_{}_iter_{}'.format(G2.graph['name'], nodes_to_delete,
                                                                         edge_frac_To_modify, i)
                    V1 = list(G1.nodes)
                    V2 = list(G2.nodes)
                    E2 = len(list(G2.edges))

                    write_leda_to_file(G1, G1_name, 'output_alignment/graphs')

                    true_A = pd.Series(data=V1, index=V1)

                    A = generate_random_alignment(V1, V2)

                    E1, Ea, Ea_hat = get_align_edges(G1, G2, A)
                    Ea_hat_star = Ea_hat
                    Ea_star = Ea
                    omega = (len(V1) * (len(V1) - 1)) / 2.0

                    nc = get_nc(A, true_A)
                    measures = get_measures(E1, Ea, Ea_hat, omega)
                    acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measures
                    row = (
                        G1_name, len(V1), E1, G2_name, len(V2), E2, omega, Ea_hat_star, Ea, Ea_hat, '', 0, '', 0,
                        '', 0, nc, tp, fp, tn, fn, ec, ics, s3, f1, acc, bm, mk, mcc, f01, f033, f3, f10)
                    rows = [row]

                    original_Ea = Ea
                    original_Ea_hat = Ea_hat
                    objective_function = measures[objective]
                    set_free_nodes = set(G2.nodes) - set(G1.nodes)
                    correct_nodes = nc * len(V1)
                    tries = 10000000
                    start = timeit.default_timer()
                    for x in range(tries):
                        if random.random() < 0:
                            v1 = random.choice(V1)
                            u1 = A[v1]
                            u2 = random.choice(list(set_free_nodes))

                            ea_add = update_Ea_change(v1, u1, u2, G1, G2, A)
                            ea_hat_add = update_Ea_hat(u1, u2, G2, set_free_nodes)

                            Ea += ea_add
                            Ea_hat += ea_hat_add

                            measures = get_measures(E1, Ea, Ea_hat, omega)
                            acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measures

                            if measures[objective] > objective_function:
                                change(A, v1, u2, set_free_nodes)
                                objective_function = measures[objective]

                                correct_nodes += update_nc(v1, u1, u2)
                            else:
                                Ea -= ea_add
                                Ea_hat -= ea_hat_add
                        else:
                            v1, v2 = random.sample(V1, 2)
                            u1, u2 = A[v1], A[v2]

                            ea_add = update_Ea_swap(v1, v2, u1, u2, G1, G2, A)

                            Ea += ea_add

                            measures = get_measures(E1, Ea, Ea_hat, omega)
                            acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measures

                            if measures[objective] > objective_function:
                                # swap(A, u1, u2, v1, v2)
                                A[v1], A[v2] = u2, u1
                                objective_function = measures[objective]

                                correct_nodes += update_nc(v1, u1, u2)
                                correct_nodes += update_nc(v2, u2, u1)
                            else:
                                Ea -= ea_add

                        nc = correct_nodes / len(V1)

                        # assert nc == get_nc(A, true_A)

                        row = (
                            G1_name, len(V1), E1, G2_name, len(V2), E2, omega, Ea_hat_star, Ea, Ea_hat, v1,
                            G1.degree[v1], u1,
                            G2.degree[u1], u2, G2.degree[u2], nc, tp, fp, tn, fn, ec, ics, s3, f1, acc, bm, mk,
                            mcc,
                            f01,
                            f033, f3, f10)
                        rows.append(row)

                        stop = timeit.default_timer()
                        print('\r sub iteration {}/{}, obj {} time: {}sec'.format(x, tries, objective_function,
                                                                                  stop - start),
                              end='')

                    df = pd.DataFrame(rows,
                                      columns=['G1', 'V1', 'E1', 'G2', 'V2', 'E2', 'omega', 'Ea_hat*', 'Ea',
                                               'Ea_hat', 'G1 changed node', 'Degree of changed node',
                                               'G2 old node', 'Degree of old node', 'G2 new node',
                                               'Degree of new node', 'NC', 'TP', 'FP', 'TN', 'FN', 'EC', 'ICS',
                                               'S3', 'F1', 'ACC', 'BM', 'MK', 'MCC', 'F0.1', 'F0.33', 'F3',
                                               'F10'])

                    output_path = 'output_alignment/%s.csv' % output_file_name
                    if os.path.isfile(output_path):
                        df.to_csv(output_path, mode='a', header=None)
                    else:
                        df.to_csv(output_path)
                    print()


def update_nc(v, u_old, u_new):
    res = 0
    if v == u_old:
        res -= 1
    if v == u_new:
        res += 1
    return res


def swap(A, u1, u2, v1, v2):
    A[v1], A[v2] = u2, u1


def generate_random_alignment(V1, V2):
    return pd.Series(data=random.sample(V2, len(V1)), index=V1)


def main():
    '''
    graphs to generate: Barabasi, Wattz, Erdos size 1000 nodes
    barabasi_albert: n: 1000, m: [9, 500] and ba(1000, 9) complement
    erdos: n: 1000, m: [8500, 250000, 491000]
    wattz: n: 1000, k: [18, 500, 982], p: 0.5

    delete [10, 50, 100] nodes
    full_iterations: 2
    shuffles: 2
    edge_frac_To_modify_list = [0.01, 0.05, 0.1]
    nodes_to_delete_list = [10, 50, 100]
    :return:
    '''

    # nx.gnm_random_graph(n, m) erdos
    # nx.watts_strogatz_graph()

    full_iterations = 5
    shuffles = 5
    # edge_frac_To_modify_list = [0.01, 0.05, 0.1]
    edge_frac_To_modify_list = [0.0]
    # nodes_to_delete_list = [10, 50, 100]
    nodes_to_delete_list = [100]

    G2s = []
    # G2s += get_barabasi_graphs()
    # G2s += get_erdos_graphs()
    # G2s += get_wattz_graphs()

    G2s += syeast0_graphs()

    # jobs = []
    # for G2 in G2s:
    #     experiment_name = '{}_iter_{}_shuffle_{}'.format(G2.graph['name'], full_iterations, shuffles)
    #     # args = ([G2], edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, shuffles, experiment_name)
    #     args = ([G2], edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, 'mcc', experiment_name)
    #     # p = Process(target=generate_alignments, args=args)
    #     p = Process(target=generate_alignments_for_objective, args=args)
    #     jobs.append(p)
    #     p.start()
    #
    # for p in jobs:
    #     p.join()
    generate_alignments(G2s, edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, shuffles, 'test')
    # generate_alignments_for_objective(G2s[:1], edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, 'mcc',
    #                                   'test_hill')

    pass


def generate_true_alignment(G, g_name, output_path):
    graph_true_alignment_path = os.path.join(output_path, '%s.align' % g_name)
    # if not os.path.exists(graph_true_alignment_path):
    #     os.makedirs(graph_true_alignment_path)
    pd.DataFrame(zip(list(G.nodes), list(G.nodes))).to_csv(graph_true_alignment_path, header=None, sep='\t',
                                                           index=False)
    # with open(graph_true_alignment_path, 'wb') as ouput_file:
    #     for node in list(G.nodes):
    #         ouput_file.write(node + '\t' + node + '\n')


def write_grapha_and_true_alignment(G1_add, G1_add_name, output_path):
    generate_true_alignment(G1_add, G1_add_name, output_path)
    write_leda_to_file(G1_add, G1_add_name, os.path.join(output_path, '{0}/'.format(G1_add_name)))


def generate_graphs_for_beta_experiments():
    output_path = 'roni_test_graphs_and_data/'
    nodes_to_delete = 100
    # np.arange(0.140, 0.162, 0.002)
    g2_names = set()
    for i in range(5):
        for target_dense in [0.017]:
            G2_name = 'syeast0'
            input_g2_file = 'networks/{0}/{0}.gw'.format(G2_name)
            G2 = nx.read_leda(input_g2_file)
            # show_graph_degree_histogram(G2)
            G2.graph['name'] = G2_name
            V2 = len(G2.nodes)
            g2_dens = nx.density(nx.read_leda(input_g2_file))
            print(g2_dens)
            # G2_base_name = G2_name.split('_')[0]
            k = int(V2 * (V2 - 1) / 2 * target_dense) - len(G2.edges)
            # # G2 = add_random_edges(G2, k)
            G2 = remove_random_edges(G2, -k)
            G2_base = '{0}_dens_{1:.4f}'.format(G2_name, nx.density(G2))
            g2_names.add(G2_base)
            write_leda_to_file(G2, G2_base, os.path.join(output_path, '{0}/'.format(G2_base)))
            G_name_template = '{}_node_del_{}_{}_edges_{}_iter_{}'
            print('G2 density {}'.format(nx.density(G2)))
            for modify in [0.30]:
                G1_name = G_name_template.format(G2_base, nodes_to_delete, 'add', modify, i)
                G1 = generate_modified_graph(G2, add_random_edges, G1_name, modify, nodes_to_delete, output_path)

                # modify = 0.5
                # G1_name = G_name_template.format(G2_name, nodes_to_delete, 'remove', modify, i)
                # generate_modified_graph(G2, remove_random_edges, G1_name, modify, nodes_to_delete, output_path)
                #
                # G1_name = G_name_template.format(G2_name, nodes_to_delete, 'add', modify, i)
                # generate_modified_graph(G2, add_random_edges, G1_name, modify, nodes_to_delete, output_path)
                #
                # G1_name = G_name_template.format(G2_name, nodes_to_delete, 'rewire', modify, i)
                # generate_modified_graph(G2, rewire_edges, G1_name, modify, nodes_to_delete, output_path)
    print(' '.join(g2_names))

def show_graph_degree_histogram(G):
    degree_sequence = np.array([d for n, d in G.degree()])
    bins = np.logspace(degree_sequence[0], degree_sequence[-1],10, base=2)
    plt.hist(degree_sequence, bins=bins)  # arguments are passed to np.histogram
    plt.title("Histogram log 2")
    plt.show()

def rewire_edges(G1, edges_to_modify):
    return nx.double_edge_swap(nx.Graph(G1), edges_to_modify, max_tries=10000)


def generate_modified_graph(G2, edge_modify_fn, G1_name, modify, nodes_to_delete, output_path):
    G1 = generate_G1(G2, nodes_to_delete, 0)
    V1 = len(G1.nodes)
    edges_to_modify = int(modify * V1)
    G1_remove = edge_modify_fn(nx.Graph(G1), edges_to_modify)
    density = nx.density(G1)
    print('G1 n: {}, m: {}, density {}, {}'.format(V1, len(G1_remove.edges), density, G1_name))
    write_grapha_and_true_alignment(G1_remove, G1_name, output_path)
    return G1


if __name__ == "__main__":
    # main()
    generate_graphs_for_beta_experiments()
