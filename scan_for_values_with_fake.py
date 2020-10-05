import sys
import os
def scan_for_values(outputHandle,inputAddress,coordinates):
    with open(inputAddress) as matchFile:
        counter = 0 
        for line in matchFile:
            line = line.strip()
            if line != 'NONE':
                allRes = line.split()
                tempList = coordinates+[str(counter)]
                for res in allRes:
                    res = res.split('/')
                    outputHandle.write('\t'.join(tempList+res)+'\n')
            counter += 1

if __name__ == '__main__':
    directory = sys.argv[1]
    output = sys.argv[2]
    # result_dict = {}
    onlyfiles = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and not f.endswith('_fake')]
    with open(output,'w') as outputFile:
        for file in onlyfiles:
            coordinates = file.split('_')
            scan_for_values(outputFile,os.path.join(directory,file),coordinates)
    with open(output+'_fakes','w') as outputFile:
        for file in onlyfiles:
            coordinates = file.split('_')
            scan_for_values(outputFile,os.path.join(directory,file+'_fake'),coordinates)