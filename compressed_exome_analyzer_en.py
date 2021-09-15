
import sys,gzip,pickle
import numpy as np

class LocalCluster:
    def __init__(self,het_members,hom_members,startbp,endbp,index,og_size) -> None:
        self.het_members = het_members
        self.hom_members = hom_members
        self.index = index
        self.startbp = startbp
        self.endbp = endbp
        self.og_size = og_size
        self.uncovered_size = og_size - (len(self.het_members)+2*len(self.hom_members))
    def compare(self,het_carriers,hom_carriers):
        intersection = len( (self.het_members & het_carriers) | (self.het_members & hom_carriers) | (self.hom_members & het_carriers) ) + (2*len(self.hom_members & hom_carriers))
        if intersection == 0:
            return None
        else:
            # 2 * len( (hom_carriers-self.het_members) | (self.) )
            return self.index,intersection,self.uncovered_size, len(self.het_members),len(self.hom_members), len(self.het_members & het_carriers),len(self.het_members &  hom_carriers), len(self.hom_members & het_carriers), len(self.hom_members & hom_carriers)
            # return self.index,intersection,len(self.members | carriers),len(self.members),self.uncovered_size
class ClusterCollection:
    def __init__(self) -> None:
        self.cluster_list = []
        self.id_list = []
        self.id_set = set()
        self.head = 0
    def load_fam(self,fam_file_addr):
        
        with open(fam_file_addr) as fam_file:
                for line in fam_file:
                    data = line.strip().split()
                    self.id_list.append(data[0])
                    if data[0] != 'NONE':
                        self.id_set.add(data[0])
        self.id_perm = pickle.load(open(fam_file_addr+'_p','rb'))
        return self.id_list,self.id_set,self.id_perm
    
    
    def load_clusters(self,cluster_file_addr):
        counter = 0 
        with gzip.open(cluster_file_addr,'rb') as cluster_file:
            for line in cluster_file:
                data = line.decode().strip().split()
                members = [item[:-2] for item in data[2:]]

                og_size = len(members)
                het_members = set()
                hom_members = set()
                for item in members:
                    if item in het_members:
                        het_members.remove(item)
                        hom_members.add(item)
                    else:
                        het_members.add(item)
                hom_members = hom_members & self.id_set
                het_members = het_members & self.id_set
                counter += 1 
                if len(het_members)+(2*len(hom_members)) < 20:
                    continue
                self.cluster_list.append(LocalCluster(het_members,hom_members,int(data[0]),int(data[1]),counter-1,og_size))

        self.cluster_list.sort(key=lambda x: (x.startbp,x.endbp))
        self.head = 0 
        return self.cluster_list[0].startbp
    
    def match_variant(self,het_carriers,hom_carriers,position):
        while self.head < len(self.cluster_list): 
            if position >self.cluster_list[self.head].endbp: #very lax filter to prevent checking to many clusters. every cluster will get a chance to stop the move forwards, however the exome addressed are also sorted
                self.head += 1
            else:
                break
        if self.head >= len(self.cluster_list):
            return -1 #its finished
        results = []
        for i in range(self.head,len(self.cluster_list)):
            if position < self.cluster_list[i].startbp: #another lax filter
                break
            elif position > self.cluster_list[i].endbp: #together with the second filter, this makes sure the cluster is covered here. 
                continue
            stat = self.cluster_list[i].compare(het_carriers,hom_carriers)
            if stat is not None:
                results.append(stat)
        return results

def write_to_file(output_file,cls_collection,chrom,position,counter,het_carriers,hom_carriers):
    stat_list = cls_collection.match_variant(het_carriers,hom_carriers,position)
    if stat_list == -1: #finished
        return -1
    elif len(stat_list) == 0: # no results
        return 0
    else:
        for item in stat_list:
            output_file.write((f'{chrom} {chrom}:{position} {counter} {len(het_carriers)} {len(hom_carriers)} '+' '.join([str(di) for di in item])+'\n').encode())
    return 1 #there were some results and they were printed to file
    
def load_exome_data(exome_file_addr,output_addr,cls_collection):
    variant_data = np.zeros(len(cls_collection.id_list)*2,dtype=np.dtype('U1'))
    headbp = cls_collection.cluster_list[0].startbp
    with open(exome_file_addr,'r') as exome_file:
        with gzip.open(output_addr,'wb') as output_file:
            with gzip.open(output_addr+'_fake','wb') as fake_output_file:
                status = 10
                counter = 0
                for line in exome_file:
                    counter+=1
                    td = line[:100].strip().split()
                    position = int(td[3])
                    if position < headbp:
                        continue
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
                        if cls_collection.id_list[index//2] == 'NONE':
                            break
                    else:
                        het_carriers = set()
                        hom_carriers = set()
                        for index in minor_allele_indices:
                            id  = cls_collection.id_list[index//2]
                            if id in het_carriers:
                                het_carriers.remove(id)
                                hom_carriers.add(id)
                            else:
                                het_carriers.add(id)
                        status = write_to_file(output_file,cls_collection,chrom,position,counter-1,het_carriers,hom_carriers)

                    for index in minor_allele_indices:
                        if cls_collection.id_list[cls_collection.id_perm[index//2]] == 'NONE':
                            break
                    else:
                        het_carriers = set()
                        hom_carriers = set()
                        for index in minor_allele_indices:
                            id  = cls_collection.id_list[cls_collection.id_perm[index//2]]
                            if id in het_carriers:
                                het_carriers.remove(id)
                                hom_carriers.add(id)
                            else:
                                het_carriers.add(id)
                        status = write_to_file(fake_output_file,cls_collection,chrom,position,counter-1,het_carriers,hom_carriers)
                    
                    if status == -1:
                        break
    
            

if __name__ == '__main__':
    cluster_file_addr = sys.argv[1]
    exome_file_addr = sys.argv[2]
    fam_file_addr = sys.argv[3]
    output_addr = sys.argv[4]

    cls_collection = ClusterCollection()
    fmobjs = cls_collection.load_fam(fam_file_addr)

    cls_collection.load_clusters(cluster_file_addr)

    load_exome_data(exome_file_addr,output_addr,cls_collection)
