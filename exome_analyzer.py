import os,sys
import numpy as np

INFO_FILE_NAME = 'info_file'
CLUSTER_FILE_NAME = '_abc.cls'
ID_MAPPER = '/sc/arion/projects/ipm2/for_roohy/UK_BB.id.mapper.NO_RELATEDS.12_18_2018.txt'


class Region:
    def __init__(self,startBP,endBP):
        self.IDList = []
        self.IDperm = None
        self.IDset = set()
        self.genoDictionary = {}
        self.fakeGenoDictionary = {}
        self.startBP = startBP
        self.endBP = endBP
    def load_fam(self,famFileAddr):
        with open(famFileAddr) as famFile:
            for line in famFile:
                data = line.strip().split()
                self.IDList.append(data[0])
                if data[0] != 'NONE':
                    self.IDset.add(data[0])
        self.IDperm = np.random.permutation(len(self.IDList))
    def load_genotyes(self,genoFileAddress):
        temp_genotype = np.zeros(len(self.IDList)*2,dtype=np.dtype('U1'))
        with open(genoFileAddress) as genoFile:
            for line in genoFile:
                data = line.strip().split()
                if int(data[3]) >= self.startBP and int(data[3]) <= self.endBP:
                    temp_genotype[:] = data[4:]
                    oneCount = np.sum(temp_genotype == '1')
                    twoCount = np.sum(temp_genotype == '2')
                    if oneCount == 0 or twoCount == 0:
                        continue
                    minorAllele = '1' if oneCount < twoCount else '2'
                    variant_name = data[0]+':'+data[3]
                    self.genoDictionary[variant_name] = set()
                    for i in range(temp_genotype.shape[0]):
                        if temp_genotype[i] == minorAllele:
                            if self.IDList[i//2] != 'NONE':
                                self.genoDictionary[variant_name].add(self.IDList[i//2])
                elif int(data[3]) > self.endBP:
                    break
    def load_doubletons(self,genoFileAddress):
        temp_genotype = np.zeros(len(self.IDList)*2,dtype=np.dtype('U1'))
        with open(genoFileAddress) as genoFile:
            for line in genoFile:
                data = line.strip().split()
                if int(data[3]) >= self.startBP and int(data[3]) <= self.endBP:
                    temp_genotype[:] = data[4:]
                    oneCount = np.sum(temp_genotype == '1')
                    twoCount = np.sum(temp_genotype == '2')
                    if oneCount == 0 or twoCount == 0:
                        continue
                    minorAllele = '1' if oneCount < twoCount else '2'
                    minorAlleleCount = oneCount if oneCount < twoCount else twoCount
                    # if minorAlleleCount!= 2 and minorAlleleCount != 3:
                    if minorAlleleCount < 10 or minorAlleleCount > 100:
                        continue
                    variant_name = data[0]+':'+data[3]
                    
                    foundFlag = True
                    fakeFoundFlag = True
                    for i in range(temp_genotype.shape[0]):
                        if temp_genotype[i] == minorAllele:
                            if self.IDList[i//2] == 'NONE':
                                foundFlag = False
                                break
                    for i in range(temp_genotype.shape[0]):
                        if temp_genotype[i] == minorAllele:
                            if self.IDList[self.IDperm[i//2]] == 'NONE':
                                fakeFoundFlag = False
                                break
                    
                    if foundFlag:
                        
                        self.genoDictionary[variant_name] = set()
                        for i in range(temp_genotype.shape[0]):
                            if temp_genotype[i] == minorAllele:    
                                self.genoDictionary[variant_name].add(self.IDList[i//2])
                    if fakeFoundFlag:
                        self.fakeGenoDictionary[variant_name] = set()
                        for i in range(temp_genotype.shape[0]):
                            if temp_genotype[i] == minorAllele:
                                self.fakeGenoDictionary[variant_name].add(self.IDList[self.IDperm[i//2]])
                elif int(data[3]) > self.endBP:
                    break
    def get_fake_ID(self,count):
        size = len(self.IDList)
        fakeSet = set()
        while count > 0:
            chosenIndex = np.random.choice(size)
            if self.IDList[chosenIndex] != 'NONE':
                fakeSet.add(self.IDList[chosenIndex])
                count -= 1
        return fakeSet
    def compare_families(self,clusterFileAddr,outputAddress):
        with open(clusterFileAddr) as clusterFile:
            with open(outputAddress,'w') as outputFile:
                with open(outputAddress+'_fake','w') as fakeOutput:
                    for line in clusterFile:
                        flag = False
                        fakeFlag = False
                        data = line.strip().split()
                        if len(data) < 2:
                            continue
                        temporary_set = set([item[:-2] for item in data])
                        temporary_intersection = len(temporary_set&self.IDset)
                        
                        for MAI in self.genoDictionary:
                            #union = len(temporary_set | self.genoDictionary[MAI])
                            #fakeUnion = len(temporary_set | self.fakeGenoDictionary[MAI])
                            intersection = len(temporary_set & self.genoDictionary[MAI])
                            
                            if intersection != 0:
                                union = len(temporary_set|self.genoDictionary[MAI])
                                outputFile.write('/'.join([str(item) for item in [MAI,len(temporary_set),
                                            len(self.genoDictionary[MAI]),temporary_intersection,
                                            intersection,union]])+'\t')
                                flag = True
                        if not flag:
                            outputFile.write('NONE')
                        
                        outputFile.write('\n')
                        
                        for MAI in self.fakeGenoDictionary:
                            fakeIntersection = len(temporary_set & self.fakeGenoDictionary[MAI])
                            
                            if fakeIntersection != 0:
                                
                                fakeUnion  = len(temporary_set|self.fakeGenoDictionary[MAI])
                                fakeOutput.write('/'.join([str(item) for item in [MAI,len(temporary_set),
                                        len(self.fakeGenoDictionary[MAI]),temporary_intersection,
                                        fakeIntersection,fakeUnion]])+'\t')
                                fakeFlag = True
                        if not fakeFlag:
                            fakeOutput.write('NONE')
                        fakeOutput.write('\n')
                        del temporary_set
    

if __name__ == '__main__':
    clusterFileAddr = sys.argv[1]
    exomeFileAddr = sys.argv[2]
    famFileAddr = sys.argv[3]
    startBP = int(sys.argv[4])
    endBP = int(sys.argv[5])
    outputAddr = sys.argv[6]
    region = Region(startBP,endBP)
    region.load_fam(famFileAddr)
    #region.load_genotyes(exomeFileAddr)
    region.load_doubletons(exomeFileAddr)
    region.compare_families(clusterFileAddr,outputAddr)

    

        