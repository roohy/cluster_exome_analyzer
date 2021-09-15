import sys,gzip,os

def scan_for_values(outputHandle,inputAddress,coordinates):
    with gzip.open(inputAddress,'rb') as matchFile:
        for line in matchFile:
            line = line.decode().strip()
            outputHandle.write((line+' '+coordinates[1]+'\n').encode())

if __name__ == '__main__':
    directory = sys.argv[1]
    output = sys.argv[2]
    # result_dict = {}
    onlyfiles = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and not f.endswith('_fake')]
    with gzip.open(output,'wb') as outputFile:
        for file in onlyfiles:
            coordinates = file.split('_')
            scan_for_values(outputFile,os.path.join(directory,file),coordinates)
    with gzip.open(output+'_fakes','wb') as outputFile:
        for file in onlyfiles:
            coordinates = file.split('_')
            scan_for_values(outputFile,os.path.join(directory,file+'_fake'),coordinates)