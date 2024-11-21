# -*- coding: utf-8 -*-
import os
import subprocess
import glob
import shutil
import argparse


parser = argparse.ArgumentParser(description="Phylogenetic analysis pipeline")
parser.add_argument("input_dir", help="Path to the input directory containing .fa files")
args = parser.parse_args()


input_dir = args.input_dir
os.makedirs('step1_mafft', exist_ok=True)

# Step 1: Sequence alignment using MAFFT
os.chdir(input_dir)
for filename in glob.glob('*.fa'):
    output_file = f"../step1_mafft/{filename}.fas"
    with open(output_file, 'w') as outfile:
        subprocess.run(['mafft', '--localpair', '--maxiterate', '1000', '--thread', '20', filename], stdout=outfile)
os.chdir('..')

# Step 2: Delete characters 'N' from sequences
os.chdir('step1_mafft')
with open('genelist.txt', 'w') as genelist:
    for filename in glob.glob('*.fas'):
        genelist.write(f"{filename}\n")

subprocess.run(['python', '/mnt/linxh/Adenophora/pipline/n_replacing_v2.py', '-g', 'genelist.txt'])
os.chdir('..')
os.makedirs('step2_delete_n', exist_ok=True)
for file in glob.glob('step1_mafft/*.new'):
    os.rename(file, f"step2_delete_n/{os.path.basename(file).replace('.new', '.fas')}")

# Step 3: TrimAl for sequence trimming
os.makedirs('step3_trimal', exist_ok=True)
for filename in glob.glob('step2_delete_n/*.fas'):
    output_file = f"step3_trimal/{os.path.basename(filename).replace('.fas', '.fasta')}"
    subprocess.run(['trimal', '-in', filename, '-out', output_file, '-gt', '0.8', '-st', '0.001'])

# Step 4: Spruceup processing
result = subprocess.run(['AMAS.py', 'concat', '-i'] + glob.glob('step3_trimal/*.fasta') + ['-f', 'fasta', '-d', 'dna'], capture_output=True, text=True)
print(result.stdout)
print(result.stderr)
if os.path.exists('concatenated.out'):
    os.rename('concatenated.out', 'concatenated.fasta')
os.makedirs('step4_1_amas', exist_ok=True)
if os.path.exists('concatenated.fasta'):
    os.rename('concatenated.fasta', 'step4_1_amas/concatenated.fasta')
if os.path.exists('partitions.txt'):
    os.rename('partitions.txt', 'step4_1_amas/partitions.txt')

os.makedirs('step4_2_spruceup', exist_ok=True)
shutil.copy('step4_1_amas/concatenated.fasta', 'step4_2_spruceup/')
shutil.copy('config.conf', 'step4_2_spruceup/')
os.chdir('step4_2_spruceup')
subprocess.run(['conda', 'run', '--no-capture-output', '--name', 'spruceup', 'python', '-m', 'spruceup', 'config.conf'])
os.chdir('..')

os.makedirs('step4_3_amas', exist_ok=True)
if os.path.exists('step4_2_spruceup/0.95_lognorms-cutoff-mo_trimmed.fasta'):
    os.rename('step4_2_spruceup/0.95_lognorms-cutoff-mo_trimmed.fasta', 'step4_3_amas/0.95_lognorms-cutoff-mo_trimmed.fasta')
if os.path.exists('step4_1_amas/partitions.txt'):
    os.rename('step4_1_amas/partitions.txt', 'step4_3_amas/partitions.txt')
os.chdir('step4_3_amas')
subprocess.run(['AMAS.py', 'split', '-f', 'fasta', '-d', 'dna', '-i', '0.95_lognorms-cutoff-mo_trimmed.fasta', '-l', 'partitions.txt', '-u', 'fasta', '-j'])
os.chdir('..')

# Step 5: TrimAl again
os.makedirs('step5_trimal', exist_ok=True)
for file in glob.glob('step4_3_amas/*.fas'):
    output_file = f"step5_trimal/{os.path.basename(file).replace('.fas', '.trimmed.fasta')}"
    subprocess.run(['trimal', '-in', file, '-out', output_file, '-gt', '0.8', '-st', '0.001'])

# Move trimmed files back
for file in glob.glob('step5_trimal/0.95_lognorms-cutoff-mo_trimmed.trimmed.fasta'):
    os.rename(file, f"step4_3_amas/{os.path.basename(file)}")

# Step 6: Exclude short sequences
os.makedirs('step6_pxs2fa', exist_ok=True)
os.chdir('step6_pxs2fa')
for file in glob.glob('../step5_trimal/*.trimmed.fasta'):
    output_file = f"{os.path.basename(file)}.fas"
    subprocess.run(['pxs2fa', '-s', file, '-o', output_file])

with open('genelist', 'w') as genelist:
    for file in glob.glob('*.fas'):
        genelist.write(f"{file}\n")

subprocess.run(['python', '/mnt/linxh/Adenophora/pipline/exclude_short_sequences.py', '-g', 'genelist', '-len', '150'])
os.chdir('..')
os.makedirs('step6_exclude', exist_ok=True)
for file in glob.glob('step6_pxs2fa/*retained.fasta'):
    os.rename(file, f"step6_exclude/{os.path.basename(file).replace('-out.trimmed.fasta.fas.retained.fasta', '.fasta')}")

# Step 7: RAxML phylogenetic analysis
os.chdir('step6_exclude')
with open('namelist.txt', 'w') as namelist:
    for file in glob.glob('*.fasta'):
        namelist.write(f"{os.path.basename(file).replace('.fasta', '')}\n")

with open('namelist.txt') as namelist:
    for name in namelist:
        name = name.strip()
        subprocess.run(['raxmlHPC-PTHREADS-SSE3', '-T', '10', '-f', 'a', '-N', '200', '-p', '12345', '-x', '12345', '-s', f"{name}.fasta", '-n', f"{name}.tre", '-m', 'GTRGAMMA'])

os.makedirs('../step7_raxml', exist_ok=True)
for file in glob.glob('RAxML_*'):
    shutil.copy(file, '../step7_raxml/')

shutil.rmtree("../step6_pxs2fa")

# Step 8: TreeShrink
with open('namelist.txt') as namelist:
    for name in namelist:
        name = name.strip()
        os.makedirs(name, exist_ok=True)
        shutil.copy(f"{name}.fasta", name)
        shutil.copy(f"RAxML_bipartitions.{name}.tre", name)
        os.chdir(name)
        os.rename(f'{name}.fasta', 'input.fasta')
        os.rename(f'RAxML_bipartitions.{name}.tre', 'input.tre')
        os.chdir('..')

os.makedirs('../step8_treeshrink_input', exist_ok=True)
for name in open('namelist.txt'):
    os.rename(name.strip(), f'../step8_treeshrink_input/{name.strip()}')

os.chdir('..')
subprocess.run(['/mnt/linxh/software/miniconda3/envs/ts/bin/run_treeshrink.py', '-i', 'step8_treeshrink_input', '-t', 'input.tre', '-a', 'input.fasta', '-o', 'step8_treeshrink_output', '-O', 'output'])

shutil.copy('step6_exclude/namelist.txt', 'step8_treeshrink_output')
os.chdir('step8_treeshrink_output')
with open('namelist.txt') as namelist:
    for name in namelist:
        name = name.strip()
        os.chdir(name)
        os.rename('output.fasta', f'{name}.fasta')
        os.rename('output.tre', f'{name}.tre')
        os.chdir('..')

os.makedirs('sequences', exist_ok=True)
os.makedirs('trees', exist_ok=True)
for fasta_file in glob.glob('./**/*.fasta', recursive=True):
    shutil.copy(fasta_file, './sequences')
for tre_file in glob.glob('./**/*.tre', recursive=True):
    shutil.copy(tre_file, './trees')

for file_path in glob.glob("../step6_exclude/RAxML*", recursive=True):
    os.remove(file_path)
