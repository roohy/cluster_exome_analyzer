import sys
from os import listdir
from os.path import isfile, join

ID_MAPPER = '/sc/orga/projects/ipm2/for_roohy/UK_BB.id.mapper.NO_RELATEDS.12_18_2018.txt'

if __name__ == '__main__':
    directory = sys.argv[1]
    onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f)) and f.endswith('.tfam')]
    idDict = {}
    with open(ID_MAPPER,'r') as mapper:
        for line in mapper:
            data = line.strip().split()
            idDict[data[2]] = data[1]
    for fileName in onlyfiles:
        print(fileName)
        with open(join(directory,fileName),'r') as famInput:
            with open(join(directory,fileName[:-5]+'.IDlist'),'w') as famOut:
                for line in famInput:
                    data = line.strip().split()
                    if data[1] in idDict:
                        famOut.write(idDict[data[1]]+'\n')
                    else:
                        famOut.write('NONE\n')
                        print("NOT FOUND "+data[1])


