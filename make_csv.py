import sys
import os
import numpy as np
import pandas as pd

def add_info_to_dataframe(array, read_length, mutation, error, kmer, assemby_data_path):
    assembly_df = pd.read_csv(assembly_data_path, sep="\t")['Contigs']
    new_array_line = [read_length, mutation, error, kmer] + list(assebly_df)[:-1]
    array.append(new_array_line)
    return array

def create_csv(directory_path, out_path):
    columns = ['read_lenght', 'mutation_rate', 'error_rate', 'kmer', 'contigs_greate_0', 'contigs_greater_1000', 'total_lenght4contigs_greater_0', 'total_length4contigs_greater_1000', 'contigs', 'largest_contig_length', 'total_length', 'GC', 'N50', 'N75', 'L50', 'L75']
    array = []
    read_lenght_paths = [os.path.join(directory_path, f) for f in os.listdir(directory_path)]
    for r in read_length_paths:
        read_lenght = int(os.path.basename(r).split('=')[1])
        mutation_rate_paths = [os.path.join(r, f) for f in os.listdir(r)]
        for m in mutation_rate_paths:
            mutation_rate = int(os.path.basename(m).split('per')[0])
            error_rate_paths = [os.path.join(m, f) for f in os.listdir(m)]
            for e in error_rate_paths:
                error_rate = int(os.path.basename(e).split('error_per')[0])
                local_velvet_test_path = os.path.join(e, 'velvet_test')
                kmer_paths = [os.path.join(local_velvet_path, f) for f in os.listdir(local_velvet_path)]
                for k in kmer_paths:
                    kmer = int(os.path.basename(k))
                    path_to_tsv = os.path.join(k, 'report.tsv')

                    add_info_to_dataframe(array, r, m, e, k, path_to_tsv)\

    my_df = pd.DataFrame(array, columns=columns)
    my_df.to_csv(out_path)

if __name__ == '__main__':
    params = sys.args[1:]
    create_csv(params[0], params[1])
