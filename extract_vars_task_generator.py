import sys,os
import subprocess

PREFIX = 'exome_'
PYTHON_SCRIPT_ADDRESS = '/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/200k/variant_extract/extract_vars.py'


if __name__ == '__main__':
    
    exome_dir = sys.argv[1]
    start_chr = int(sys.argv[2])
    end_chr = int(sys.argv[3])
    sub_dir = sys.argv[4]
    output_addr = sys.argv[5]
    for i in range(start_chr,end_chr+1):    
        subfile_addr = os.path.join(sub_dir,f'subfile_chr{i}')
        with open(subfile_addr,'w') as subfile:
            subfile.write(f'#BSUB -J exome_chr{i}\n')
            subfile.write('#BSUB -P acc_ipm2\n#BSUB -q premium\n#BSUB -n 1\n'+
                     '#BSUB -R "span[hosts=1] affinity[core(2, same=socket, exclusive=(socket, injob))]"\n'+
                     '#BSUB -R rusage[mem=3000]\n#BSUB -W 24:00\n')
            subfile.write('#BSUB -o '+os.path.join(sub_dir,f'chr{i}_job.stdout')+'\n')
            subfile.write('#BSUB -e '+os.path.join(sub_dir,f'chr{i}_job.stderr')+'\n')
            subfile.write('#BSUB -L /bin/bash\n')
            subfile.write('module load anaconda3\n')
            args = ['python', PYTHON_SCRIPT_ADDRESS,
                    os.path.join(exome_dir,'output_chr'+str(i)+'.tped'),
                    os.path.join(exome_dir,'output_chr'+str(i)+'.id_list'),
                    os.path.join(output_addr,f'chr_{i}_vars')]
            subfile.write(' '.join(args))
        subprocess.run('bsub < '+subfile_addr,shell=True)
        print(subfile_addr)


