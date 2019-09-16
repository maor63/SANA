import os

import pandas as pd

base_dir = 'y0_y2_0.001_0.015_dens_s3_ec_ics_fbetahash_5min_25_iter'
run_times_dirs_names = os.listdir(base_dir)

measures = ['ec', 'ics', 's3', 'lccs', 'sec', 'nc', 'wacc', 'mcc', 'bm', 'mk', 'fbetastar', 'fbetahash']
header = ['g1', 'g1_n', 'g1_edges', 'density_g1', 'g2', 'g2_n', 'g2_edges', 'density_g2', 'method',
          'run_time'] + measures
csv_output_rows = []
for run_time in run_times_dirs_names:
    objective_functions = os.listdir(base_dir + '/' + run_time)
    for objective_function in objective_functions:
        experiments = os.listdir(base_dir + '/' + run_time + '/' + objective_function)
        for experiment in experiments:
            try:
                output_file = open(
                    base_dir + '/' + run_time + '/' + objective_function + '/' + experiment + '/sana_' + objective_function.lower() + '.out')
                res_row = []
                score_results = {}
                for row in output_file:
                    if row.startswith('G1:'):
                        g1_name = row.split(' ')[1].replace('\n', '')
                        g1_nodes = output_file.next().split(' ')[-1].replace('\n', '')
                        g1_edges = output_file.next().split(' ')[-1].replace('\n', '')
                        nodes = int(g1_nodes)
                        density_g1 = round(float(g1_edges) / (nodes * (nodes - 1) / 2), 4)
                        res_row += [g1_name, g1_nodes, g1_edges, density_g1]
                    elif row.startswith('G2:'):
                        g2_name = row.split(' ')[1]
                        g2_nodes = output_file.next().split(' ')[-1].replace('\n', '')
                        g2_edges = output_file.next().split(' ')[-1].replace('\n', '')
                        nodes = int(g2_nodes)
                        density_g2 = round(float(g2_edges) / (nodes * (nodes - 1) / 2), 4)
                        res_row += [g2_name, g2_nodes, g2_edges, density_g2]
                    elif row.startswith('Execution time:'):
                        execution_time = row.split(' ')[-1].replace('\n', '')
                        res_row += [execution_time]
                    elif row.startswith('Method:'):
                        method_name = row.split(' ')[1].replace('\n', '')
                        res_row += [objective_function]
                    elif row.startswith('Scores:'):
                        for score_row in output_file:
                            if score_row == ' \n' or score_row == '\n':
                                break
                            score_name, score_result = score_row.split(' ')
                            score_results[score_name.replace(':', '')] = score_result.replace('\n', '')
                for measure in measures:
                    res_row.append(score_results[measure])
                csv_output_rows.append(res_row)
            except Exception as e:
                print(e)

df = pd.DataFrame(data=csv_output_rows, columns=header)
df.to_csv('{}.csv'.format(base_dir))
