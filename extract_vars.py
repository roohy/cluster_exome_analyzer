
import sys,gzip,pickle
import numpy as np

class FamObject:
    def __init__(self) -> None:
        
        self.id_list = []
        self.id_set = set()
        
    def load_fam(self,fam_file_addr):
        
        with open(fam_file_addr) as fam_file:
                for line in fam_file:
                    data = line.strip().split()
                    self.id_list.append(data[0])
                    if data[0] != 'NONE':
                        self.id_set.add(data[0])
        self.id_perm = np.random.permutation(len(self.id_list))

        pickle.dump(self.id_perm,open(fam_file_addr+'_p','wb'))
    
    
def write_to_file(output_file,chrom,position,counter,het_carriers,hom_carriers):
    
    output_file.write((f'{chrom} {chrom}:{position} {counter} {len(het_carriers)} {len(hom_carriers)}\n').encode())
    
    
def load_exome_data(exome_file_addr,output_addr,famobj):
    variant_data = np.zeros(len(famobj.id_list)*2,dtype=np.dtype('U1'))
    
    with open(exome_file_addr,'r') as exome_file:
        with gzip.open(output_addr,'wb') as output_file:
            status = 10
            counter = 0 
            for line in exome_file:
                counter += 1
                td = line[:100].strip().split()
                position = int(td[3])
                data = line.strip().split()
                # position = int(data[3])
                chrom = int(data[0])
                variant_data[:] = data[4:]
                one_count = np.sum(variant_data == '1')
                two_count = np.sum(variant_data == '2')
                if one_count == 0 or two_count == 0:
                    continue
                minor_allele = '1' if one_count < two_count else '2'
                MAC = one_count if one_count < two_count else two_count
                if MAC > 4000 or MAC < 10:
                    continue
                # variant_name = data[0]+':'+data[3]
                # found_all = True
                # fake_found_all = True
                minor_allele_indices = np.where(variant_data == minor_allele)[0]
                for index in minor_allele_indices:
                    if famobj.id_list[index//2] == 'NONE':
                        break
                else:
                    het_carriers = set()
                    hom_carriers = set()
                    for index in minor_allele_indices:
                        id  = famobj.id_list[index//2]
                        if id in het_carriers:
                            het_carriers.remove(id)
                            hom_carriers.add(id)
                        else:
                            het_carriers.add(id)
                    write_to_file(output_file,chrom,position,counter-1,het_carriers,hom_carriers)

                    
                    
            

if __name__ == '__main__':
    exome_file_addr = sys.argv[1]
    fam_file_addr = sys.argv[2]
    output_addr = sys.argv[3]

    famobj = FamObject()
    famobj.load_fam(fam_file_addr)

    load_exome_data(exome_file_addr,output_addr,famobj)
