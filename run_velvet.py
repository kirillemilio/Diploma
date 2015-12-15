import sys
import os
import os.path.join
import subprocess
from multiprocessing import Process

def mult_velvet_run(param_list):
    velvet_pr = []
    for item in param_list:
        pr = Process(target=run_velvet, args=item)
        pr.start()
        velvet_pr.append(pr)
    return velvet_pr

def run_velvet(velvet_path, quast_path, reads, out_velvet_dir, out_velvet_test_dir, current_hash):
    subprocess.call([join(velvet_path, 'velveth'), out_velvet_dir, current_hash, '-fastq', reads])
    subprocess.call([join(velvet_path, 'velvetg'), out_velvet_dir])

def join_velvet_threads(threads):
    for thread in threads:
        thread.join()
    return
