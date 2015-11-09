import subprocess
import os
import sys
from os.path import join
# from os import isfile
# from os import listdir
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import random

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

make_mutations('/home/kirill/git_assemblers/bacteria.fasta', '/home/kirill/git_assemblers/bacteria_alt.fasta', 20)
