import sys,os
import subprocess

PREFIX = 'exome_'
PYTHON_SCRIPT_ADDRESS = '/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/exome_analyzer_compressed.py'


if __name__ == '__main__':
    libd_dir = sys.argv[1]
    exome_dir = sys.argv[2]
    start_chr = int(sys.argv[3])
    end_chr = int(sys.argv[4])
    sub_dir = sys.argv[5]
    output_addr = sys.argv[6]
    for i in range(start_chr,end_chr+1):
        graph_dir = os.path.join(os.path.join(libd_dir,str(i)),'vars')
        windows = [(os.path.join(graph_dir,f),int(f[:-3])) for f in os.listdir(graph_dir) if f.endswith('.gz')]
        for window in windows:
            subfile_addr = os.path.join(sub_dir,f'subfile_{i}_{window[1]}')
            with open(subfile_addr,'w') as subfile:
                subfile.write(f'#BSUB -J exome_cls_{i}_{window[1]}\n')
                subfile.write('#BSUB -P acc_ipm2\n#BSUB -q premium\n#BSUB -n 1\n'+
                         '#BSUB -R "span[hosts=1] affinity[core(2, same=socket, exclusive=(socket, injob))]"\n'+
                         '#BSUB -R rusage[mem=100000]\n#BSUB -W 10:00\n')
                subfile.write('#BSUB -o '+os.path.join(sub_dir,f'{i}_{window[1]}.stdout')+'\n')
                subfile.write('#BSUB -e '+os.path.join(sub_dir,f'{i}_{window[1]}.stderr')+'\n')
                subfile.write('#BSUB -L /bin/bash\n')
                subfile.write('module load anaconda3\n')
                args = ['python', PYTHON_SCRIPT_ADDRESS,window[0],
                        os.path.join(exome_dir,'output_chr'+str(i)+'.tped'),
                        os.path.join(exome_dir,'output_chr'+str(i)+'.id_list'),
                        os.path.join(output_addr,f'{i}_{window[1]}')]
                subfile.write(' '.join(args))
            subprocess.run('bsub < '+subfile_addr,shell=True)
            print(subfile_addr)


