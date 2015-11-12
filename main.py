import subprocess
import os
import sys
from os.path import join
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import random
import shutil

def unite2fastq(file1, file2, destination):
    out_file = open(destination, 'w')
    f1 = open(file1)
    f2 = open(file2)
    shutil.copyfileobj(file1, out_file)
    shutil.copyfileobj(file2, out_file)
    out_dile.close()
    file1.close()
    file2.close()

def make_mutations(fasta_path, fasta_out, mutation_rate):
    fasta = list(SeqIO.parse(fasta_path, 'fasta'))[0]
    st = list(str(fasta.seq))
    mutation_position_list = random.sample(range(len(st)), int(mutation_rate * len(st) / 100))
    for pos in mutation_position_list:
        st[pos] = possible_letters = random.choice(filter(lambda x: x != st[pos], ['A', 'T', 'G', 'C']))
    s = ''.join(st)
    print(float(len(mutation_position_list)) / len(s))
    new_seq = Seq(s, IUPAC.unambiguous_dna)
    fasta.seq = new_seq
    SeqIO.write(fasta, fasta_out, 'fasta')

def run_dwg(dwg_path, input_fasta, out_dir, prefix, error_rate, read_length, mean_coverage=30, random_read_prob=0, mutation_rate=0, distance_between_read_pairs=1000):
    subpreocess.call([dwg_path, '-e', str(error_rate), '-E', str(error_rate), '-d', str(distance_between_read_pairs), '-C', str(mean_coverage), '-1', str(read_length), '-2', str(read_length), '-y', str(random_read_prob), '-r', str(mutation_rate), input_fasta, join(out_dir, prefix)])
    return join(out_dir, prefix + '.bfast.fastq')

def generate_dirs(input_genome_fasta, out_dir, read_length_list, mutations_list, error_list):
    if(os.path.exists(out_dir)):
        shutil.rmtree(out_dir)
    os.mkdir(out_dir)
    for read_length in read_length_list:
        read_length_prefix = join(out_dir, 'READ_LENGTH=' + str(read_length))
        os.mkdir(read_length_prefix)
        for mutation in mutations_list:
            mutation_prefix = join(read_length_prefix, str(mutation) + "per100")
            os.mkdir(mutation_prefix)
            for error in error_list:
                error_prefix = join(mutation_prefix, str(error) + "error_per_100")
                os.mkdir(error_prefix)
                shutil.copy(input_genome_fasta, os.path. os.path.basename(input_genome_fasta))



# make_mutations('/home/kirill/git_assemblers/bacteria.fasta', '/home/kirill/git_assemblers/bacteria_alt.fasta', 20)
