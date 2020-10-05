import sys
import os


if __name__ == '__main__':
    directory = sys.argv[1]
    output = sys.argv[2]
    # result_dict = {}
    onlyfiles = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and not f.endswith('_fake')]
    with open(output,'w') as outputFile:
        # resultArray = [[],[],[]]
        for file in onlyfiles:
    
            coordinates = file.split('_')
            # if coordinates[0] not in result_dict:
            # result_dict[coordinates[0]] = {coordinates[1]:{}}
            with open(os.path.join(directory,file)) as matchFile:
                counter = 0
                for line in matchFile:
                    line = line.strip()
                    if line != 'NONE':
                        allRes = line.split()
                        temp_list = coordinates+[str(counter)]
                        for res in allRes:
                            res = res.split('/')
                            outputFile.write('\t'.join(temp_list+res) + '\n')
                    counter += 1