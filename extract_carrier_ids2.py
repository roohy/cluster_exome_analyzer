import sys,pickle
import numpy as np 
import pandas as pd
import pyarrow.feather as feather
import extract_family2 as ef
if __name__ == '__main__':
    df = feather.read_feather('/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/200k/0.75rec_0.1fp_30us_20c_25_10_4000_en.feather')
    

    cov = ef.Covariates()
    cov.add_from_file('/sc/arion/projects/ipm2/for_roohy/UKBB.age_sex.pcs.sorted.covar')
    # cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/hdl_chol.ukbb.mapped.pheno')
    
    cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/urae.ukbb_mapped.pheno')
    # cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/standing_height.ukbb.mapped.pheno')
    context = ef.ExomeContext('/sc/arion/projects/ipm2/roohy/ukbb_ibd/global_network/inon/output_chr1.id_list')
    fam_ext = ef.FamilyExtractor('/sc/arion/projects/ipm2/roohy/ukbb_ibd/umapping_25_i2_cls/')
    car_ext = ef.CarrierExtractor('/sc/arion/projects/ipm2/roohy/ukbb_ibd/global_network/inon/',context)
    
    #FIRST TIME
    #carriers = car_ext.get_variables(df['var_name'],df['var_index'])

    #pickle.dump(carriers,open('/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/200k/carriers_75rec1fp30_20_4k.pkl','wb'))

    #second time
    carriers = pickle.load(open('/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/200k/carriers_75rec1fp30_20_4k.pkl','rb'))
    
    # filtered_carriers = ef.filter_results(tops,carriers,2,3)
    glm = ef.SingleGLM()
    glm.set_cov(cov.pandas,context)
    pvs = ef.test_pipe(carriers,df,fam_ext,glm)
    pickle.dump(pvs,open('pvs_75rec10fp20_20_urae.pkl','wb'))

    cov = ef.Covariates()
    cov.add_from_file('/sc/arion/projects/ipm2/for_roohy/UKBB.age_sex.pcs.sorted.covar')
    cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/creatinine.ukbb.mapped.pheno')
    # cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/chol.ukbb.mapped.pheno')

    glm = ef.SingleGLM()
    glm.set_cov(cov.pandas,context)
    pvs = ef.test_pipe(carriers,df,fam_ext,glm)
    pickle.dump(pvs,open('pvs_75rec10fp20_20_creatinine.pkl','wb'))