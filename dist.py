import argparse
import re
import numpy as np

parser = argparse.ArgumentParser('python dist.py')
parser.add_argument('depth', metavar='N', type=float, nargs='+',  help='A set of depths for which to sample')
parser.add_argument('--input', '-i', type=str, help='An input fastq file')
parser.add_argument('--output', '-o', type=str, help='An output directory for subsampled files')

args = parser.parse_args()

depths = args.depth

file = args.input
output = re.sub(r'/$', '', args.output)+r'/'

for sample in depths:
    if '2P' in file:
        continue
    else:
        with open(file, 'r') as f:
            r1 = f.readlines()
        with open(re.sub(r'1P.fastq', r'2P.fastq', file), 'r') as f:
            r2 = f.readlines()
        outfileName = output+file.split('/')[-1]+'_'+str(sample)
        with open(outfileName, 'w') as fd1, open(re.sub(r'1P.fastq', r'2P.fastq', outfileName), 'w') as fd2:
            total = len(r1) / 4
            idx = np.random.choice(total, total*sample, replace=False)
            
            print(outfileName)
            print(len(r1), total, total*sample)
            print(len(r2), total, total*sample)
            if len(r1) != len(r2):
                print('Warning: Unequal file lengths detected for '+outfileName+', check values')

            for i in np.nditer(idx):
                fd1.write(r1[i*4])
                fd1.write(r1[i*4 + 1])
                fd1.write(r1[i*4 + 2])
                fd1.write(r1[i*4 + 3])
                fd2.write(r2[i*4])
                fd2.write(r2[i*4 + 1])
                fd2.write(r2[i*4 + 2])
                fd2.write(r2[i*4 + 3])
