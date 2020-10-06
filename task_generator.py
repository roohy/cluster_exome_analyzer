import sys,os
import subprocess

PREFIX = 'exome_'
PYTHON_SCRIPT_ADDRESS = '/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/exome_analyzer.py'


if __name__ == '__main__':
    directory = sys.argv[1]
    exomeDirectory = sys.argv[2]
    startChr = int(sys.argv[3])
    endChr = int(sys.argv[4])
    subDir = sys.argv[5]
    outputAddr = sys.argv[6]
    for i in range(startChr,endChr+1):
        basePath = os.path.join(directory,str(i),'graphs')
        with open(os.path.join(basePath,'info_file')) as infoFile:
            for line in infoFile:
                data = line.strip().split()
                #data = [int(item) for item in data]
                subfileAddr = os.path.join(subDir,'subfile_'+str(i)+'_'+data[0])
                with open(subfileAddr,'w') as subfile:
                    subfile.write('#BSUB -J exome_cls_'+str(i)+'_'+data[0]+'\n')
                    subfile.write('#BSUB -P acc_ipm2\n#BSUB -q premium\n#BSUB -n 1\n'+
                         '#BSUB -R "span[hosts=1] affinity[core(2, same=socket, exclusive=(socket, injob))]"\n'+
                         '#BSUB -R rusage[mem=100000]\n#BSUB -W 10:00\n')
                    subfile.write('#BSUB -o '+os.path.join(subDir,str(i)+'_'+data[0]+'.stdout')+'\n')
                    subfile.write('#BSUB -e '+os.path.join(subDir,str(i)+'_'+data[0]+'.stderr')+'\n')
                    subfile.write('#BSUB -L /bin/bash\n')
                    subfile.write('module load anaconda3\n')
                    args = ['python', PYTHON_SCRIPT_ADDRESS,os.path.join(basePath,data[0],'_abc.cls'),
                        os.path.join(exomeDirectory,'exome_'+str(i)+'.tped'),
                        os.path.join(exomeDirectory,'exome_'+str(i)+'.IDlist'),
                        data[1],data[2],os.path.join(outputAddr,str(i)+'_'+data[0])]
                    subfile.write(' '.join(args))
                subprocess.run('bsub < '+subfileAddr,shell=True)
                print(subfileAddr)


