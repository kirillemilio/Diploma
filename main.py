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
# import subprocess
from multiprocessing import Process

def unite2fastq(file1, file2, destination):
    out_file = open(destination, 'w')
    f1 = open(file1)
    f2 = open(file2)
    shutil.copyfileobj(f1, out_file)
    shutil.copyfileobj(f2, out_file)
    out_file.close()
    f1.close()
    f2.close()

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
    subprocess.call([dwg_path, '-e', str(error_rate), '-E', str(error_rate), '-d', str(distance_between_read_pairs), '-C', str(mean_coverage), '-1', str(read_length), '-2', str(read_length), '-y', str(random_read_prob), '-r', str(mutation_rate), input_fasta, join(out_dir, prefix)])
    return join(out_dir, prefix + '.bfast.fastq')

def mult_velvet_run(param_list, threads_limit=100):
    velvet_pr = []
    for item in param_list:
        pr = Process(target=run_velvet, args=item)
        pr.start()
        velvet_pr.append(pr)
        if(len(velvet_pr) >= threads_limit):
            join_velvet_threads(velvet_pr)
            velvet_pr = []
    return velvet_pr

def run_velvet(velvet_path, quast_path, reads, out_velvet_dir, out_velvet_test_dir, current_hash):
    os.mkdir(out_velvet_dir)
    os.mkdir(out_velvet_test_dir)
    subprocess.Popen([join(velvet_path, 'velveth'), out_velvet_dir, str(current_hash), '-fastq', reads], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.Popen([join(velvet_path, 'velvetg'), '-min_contig_lght', str(10), out_velvet_dir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.Popen([join(quast_path, 'quast.py'), '-o', out_velvet_test_dir, join(out_velvet_dir, 'contigs.fa')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("running quast")
    shutil.rmtree(out_velvet_dir)
    compress_quast_results(out_velvet_test_dir)

def join_velvet_threads(threads):
    for thread in threads:
        thread.join()
    return

def compress_quast_results(out_velvet_test_dir):
    pass

def generate_dirs(input_genome_fasta='/home/kemelyanov/Diploma/Initial_information/bacteria.fasta', out_dir='/home/kemelyanov/Diploma/Out', dwg_path='/home/kemelyanov/Diploma/Program/DWGSIM/dwgsim', v_path='/home/kemelyanov/Diploma/Program/velvet', s_path='', q_path='/home/kemelyanov/Diploma/Program/Quast', read_length_list=[400], mutations_list=[0, 5, 10], error_list=[0]):
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

                original_fasta = join(error_prefix, os.path.basename(input_genome_fasta))
                mutated_fasta = join(error_prefix, '_mutated.'.join(os.path.basename(input_genome_fasta).split('.')))
                shutil.copy(input_genome_fasta, original_fasta)
                make_mutations(input_genome_fasta, mutated_fasta, mutation)

                dwg_original_temp = join(error_prefix, 'temp1')
                dwg_mutated_temp = join(error_prefix, 'temp2')
                os.mkdir(dwg_original_temp)
                os.mkdir(dwg_mutated_temp)

                run_dwg(dwg_path, original_fasta, dwg_original_temp, 'temp', error, read_length)
                run_dwg(dwg_path, mutated_fasta, dwg_mutated_temp, 'temp', error, read_length)

                united_reads_path = join(error_prefix, 'reads.fastq')
                unite2fastq(join(dwg_original_temp, 'temp.bfast.fastq'), join(dwg_mutated_temp, 'temp.bfast.fastq'), united_reads_path)

                shutil.rmtree(dwg_original_temp)
                shutil.rmtree(dwg_mutated_temp)

                os.remove(original_fasta)
                os.remove(mutated_fasta)

                core_directory = error_prefix
                velvet_out = join(core_directory, 'velvet_out')
                if(not os.path.exists(velvet_out)):
                    os.mkdir(velvet_out)

                velvet_test_out = join(core_directory, 'velvet_test')
                if(not os.path.exists(velvet_test_out)):
                    os.mkdir(velvet_test_out)

                hash_lengths = [h for h in range(20, read_length + 1) if h % 2 != 0]
                velvet_paths = [join(velvet_out, str(h)) for h in hash_lengths]
                quast_paths = [join(velvet_test_out, str(h)) for h in hash_lengths]
                quast = [q_path] * len(hash_lengths)
                velvet = [v_path] * len(hash_lengths)
                reads4velvet = [united_reads_path] * len(hash_lengths)
                params = zip(velvet, quast, reads4velvet, velvet_paths, quast_paths, hash_lengths)
                # print(params)
                # input()

                threads = mult_velvet_run(params)
                join_velvet_threads(threads)

generate_dirs()
# make_mutations('/home/kirill/git_assemblers/bacteria.fasta', '/home/kirill/git_assemblers/bacteria_alt.fasta', 20)
