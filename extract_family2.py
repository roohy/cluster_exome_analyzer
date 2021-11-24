import sys,os,gzip
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import statsmodels.api as sm
import gzip
from scipy.stats import ttest_ind,chisquare,fisher_exact

TPED_PRE = 'output_chr'
TPED_SUFF = '.tped'
MAC_UPPER = 200
MAC_LOWER = 10


class SingleGLM:
    def __init__(self) -> None:
        self.covariates = None
    def set_cov(self,covariates,context):
        self.context = context
        
        self.full_covariates = covariates
        self.covariates = covariates[covariates[covariates.columns[1]].isin(context.id_set)]

    def do_glm(self,het_members,hom_members,remove_list,add_mode=False):

        covariates = self.covariates[~self.covariates[self.covariates.columns[1]].isin(remove_list)]
        if add_mode:
            addons = self.full_covariates[self.full_covariates[self.full_covariates.columns[1]].isin((hom_members|het_members)-self.context.id_set)]
            covariates = covariates.append(addons,ignore_index=True)
        covariates['target'] = 0
        covariates.loc[covariates[covariates.columns[1]].isin(het_members),'target'] = 1 
        covariates.loc[covariates[covariates.columns[1]].isin(hom_members),'target'] = 2
        x = covariates.drop([covariates.columns[0],covariates.columns[1],'pheno'],axis=1)
        '''x = covariates.drop([covariates.columns[0],covariates.columns[1],covariates.columns[2],
            covariates.columns[3],covariates.columns[5],covariates.columns[6],
            covariates.columns[7],covariates.columns[8],'target','pheno'],axis=1)'''
        model = sm.OLS(covariates['pheno'],x ).fit()
        #print(model.summary())
        result = model.pvalues['target']
        oddef = model.params['target']

        return result,oddef
    def do_two(self,het_carriers,hom_carriers,het_members,hom_members):
        het_unc = het_members - (het_carriers | hom_carriers)
        hom_unc = hom_members - (het_carriers | hom_carriers)
        remove_list = het_unc | hom_unc
        pval,oddef = self.do_glm(het_carriers,hom_carriers,remove_list)
        test_pval,test_oddef = self.do_glm(het_unc,hom_unc,het_carriers|hom_carriers,add_mode=True)
        return [[pval,oddef],[test_pval,test_oddef]]
    def return_predictions(self,addon_members):
        covariates = self.covariates
        addons = self.full_covariates[self.full_covariates[self.full_covariates.columns[1]].isin(addon_members-self.context.id_set)]
        covariates = covariates.append(addons,ignore_index=True)
        
        x = covariates.drop([covariates.columns[0],covariates.columns[1],'pheno'],axis=1)
        '''x = covariates.drop([covariates.columns[0],covariates.columns[1],covariates.columns[2],
            covariates.columns[3],covariates.columns[5],covariates.columns[6],
            covariates.columns[7],covariates.columns[8],'target','pheno'],axis=1)'''
        model = sm.OLS(covariates['pheno'],x ).fit()
        #print(model.summary())
        
        return model,covariates

class Covariates(object):

    def __init__(self) -> None:
        self.pandas  = None
        self.covar_header = None
        #self.data  = None # A list of Numpy objects 

    def add_from_file(self,file_address,separator='\s+',header=None):
        self.pandas = pd.read_csv(file_address,sep=separator,header=header)
        if header is not None:
            self.covar_header = True
        else:
            self.covar_header = False

    def add_pheno(self,path,separator='\s+',header=None):
        self.pheno = pd.read_csv(path,sep=separator,header=header, names=['fam_id','id','pheno'])
        self.pandas = self.pandas.join(self.pheno.set_index('id'),on=self.pandas.columns[1],how='inner')
        self.pandas = self.pandas.drop(['fam_id'],axis=1)
        self.pandas['ones'] = np.zeros((self.pandas.shape[0]))
        self.pandas['ones'] = 1
        self.pandas = self.pandas.dropna()
    def filter_population(self,population):
        self.pandas = self.pandas[self.pandas[1].isin(population.sample_set)]
        print(self.pandas.shape)
        #df = df[~df['date'].isin(a)]
    def filter_gender(self,gender_code):
        self.pandas = self.pandas[self.pandas[self.pandas.columns[2]]==gender_code]
        self.pandas = self.pandas.drop(self.pandas.columns[2],axis=1)
    def save_file(self,path):
        #self.pandas = (self.pandas-self.pandas.mean())/self.pandas.std()
        self.pandas['ones'] = np.zeros((self.pandas.shape[0]))
        self.pandas['ones'] = 1
        self.pandas = self.pandas.dropna()
        self.pandas.to_csv(path,sep=' ',index=False)
    def randomize_ids(self):
        self.pandas[self.pandas.columns[1]] = np.random.permutation(self.pandas[self.pandas.columns[1]].values)




class ExomeContext:
    def __init__(self,fam_file_addr) -> None:
        self.id_list = []
        self.id_set = set()
        with open(fam_file_addr) as fam_file:
                for line in fam_file:
                    data = line.strip().split()
                    self.id_list.append(data[0])
                    if data[0] != 'NONE':
                        self.id_set.add(data[0])

class FamilyExtractor:
    def __init__(self,directory) -> None:
        self.directory = directory
    
    def extract_family(self,chr_num,file_number,index):
        index = int(index)
        file_addr = os.path.join(os.path.join(self.directory,str(int(chr_num))),'vars',f'{int(file_number)}.gz')
        with gzip.open(file_addr,'rb') as cls_file:
            counter = 0 
            for line in cls_file:
                if counter != index:
                    counter += 1
                    continue
                else:
                    data = line.decode().strip().split()
                    members = [item[:-2] for item in data[2:]]
                    het_members = set()
                    hom_members = set()
                    for item in members:
                        if item in het_members:
                            het_members.remove(item)
                            hom_members.add(item)
                        else:
                            het_members.add(item)
                    return het_members,hom_members
        raise ValueError
        
def filter_results(df,carrier_list,het_ind,hom_ind):
    result = {}
    for key in carrier_list:
        result[key] = []
        for item in carrier_list[key]:
            if (len(item[0]) == df[df['var_name']==key].values[0,het_ind]) and (len(item[1]) == df[df['var_name']==key].values[0,hom_ind]):
                result[key].append(item)
    return result

    
class CarrierExtractor:
    def __init__(self,directory,context):
        self.directory = directory
        self.context = context
    def get_variables(self,var_name_list,var_index_list):
        if len(var_name_list) == 0:
            raise ValueError
        var_list = zip(var_name_list,var_index_list)
        var_list = [item[0].split(':')+[item[1]] for item in var_list]
        var_list.sort(key=lambda x:(int(x[0]),int(x[2])))
        carriers_list = {}
        variant_data = np.zeros(len(self.context.id_list)*2,dtype=np.dtype('U1'))
        # var_chr,var_position = var_name.split(':')
        current_chr = '0' #var_list[0][0]
        loc=-1
        for ind,var in enumerate(var_list):
            var_name = f'{var[0]}:{var[1]}'
            # carriers_list[var_name] =[]
            print(f'looking for var {var[1]}')
            if var[0] != current_chr:
                print('Chaging the chromosome!')
                chr_dir = os.path.join(self.directory,f'{TPED_PRE}{var[0]}{TPED_SUFF}')
                current_chr = var[0]
                chr_file = open(chr_dir,'r')  
                loc=-1
            for line in chr_file:
                loc += 1
                if loc == int(var[2]):
                    data = line.strip().split()
                    variant_data[:] = data[4:]
                    carriers = self.extract_carriers(variant_data)
                    if carriers is not None:
                        carriers_list[var_name]=carriers
                elif loc > int(var[2]):
                    break
        return carriers_list

    def extract_carriers(self,variant_data):
        one_count = np.sum(variant_data == '1')
        two_count = np.sum(variant_data == '2')
        if one_count == 0 or two_count == 0:
            return None
        minor_allele = '1' if one_count < two_count else '2'
        MAC = one_count if one_count < two_count else two_count
        if MAC > MAC_UPPER or MAC < MAC_LOWER:
            return None
        minor_allele_indices = np.where(variant_data == minor_allele)[0]    
        for index in minor_allele_indices:
            if self.context.id_list[index//2] == 'NONE':
                return None
        het_carriers = set()
        hom_carriers = set()
        for index in minor_allele_indices:
            id  = self.context.id_list[index//2]
            if id in het_carriers:
                het_carriers.remove(id)
                hom_carriers.add(id)
            else:
                het_carriers.add(id)
        return [het_carriers,hom_carriers]

def test_pipe(carriers,tops,fam_ext,glm):
    #filtered_carriers = filter_results(tops,carriers,2,3)
    # for key in filtered_carriers:
    #     #4 and 13
    #     line = tops[tops['var_name']==key].values[0,:]
    #     cls_ind = int(line[4])
    #     file_ind = int(line[13])
    #     chr_num = int(line[0])
    #     het_members,hom_members = fam_ext.extract_family(chr_num,file_ind,cls_ind)
    #     flag = False
    #     for item in filtered_carriers[key]:
    #         het_carriers = item[0]
    #         hom_carriers = item[1]
    #         intersection = len( (het_members & het_carriers) | (het_members & hom_carriers) | (hom_members & het_carriers) ) + (2*len(hom_members & hom_carriers))
    #         if flag and intersection==int(line[5]):
    #             print('dup',line)
    pvs = {}
    counter =-1 
    for key in carriers:
        counter += 1
        #4 and 13
        pvs[key] = []
        line = tops[tops['var_name']==key].values[0,:]
        cls_ind = int(line[5])
        file_ind = int(line[14])
        chr_num = int(line[0])
        het_members,hom_members = fam_ext.extract_family(chr_num,file_ind,cls_ind)
        flag = False
        
        het_carriers = carriers[key][0]
        hom_carriers = carriers[key][1]
        intersection = len( (het_members & het_carriers) | (het_members & hom_carriers) | (hom_members & het_carriers) ) + (2*len(hom_members & hom_carriers))
        if intersection != int(line[6]):
            print(f'{key} : {key} Pfffft')
            
        pvs[key].append(glm.do_two(het_carriers,hom_carriers,het_members,hom_members))
    return pvs
    
        
    
    