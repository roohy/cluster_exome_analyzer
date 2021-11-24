import sys,pickle
import numpy as np 
import pandas as pd
import pyarrow.feather as feather
import extract_family as ef
if __name__ == '__main__':
    df = feather.read_feather('/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/200k/90_rec_15_fp_10_10.feather')
    

    cov = ef.Covariates()
    cov.add_from_file('/sc/arion/projects/ipm2/for_roohy/UKBB.age_sex.pcs.sorted.covar')
    # cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/hdl_chol.ukbb.mapped.pheno')
    cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/standing_height.ukbb.mapped.pheno')
    context = ef.ExomeContext('/sc/arion/projects/ipm2/roohy/ukbb_ibd/global_network/inon/output_chr1.id_list')
    fam_ext = ef.FamilyExtractor('/sc/arion/projects/ipm2/roohy/ukbb_ibd/umapping_25_i2_cls/')
    car_ext = ef.CarrierExtractor('/sc/arion/projects/ipm2/roohy/ukbb_ibd/global_network/inon/',context)
    carriers = car_ext.get_variables(df['var_name'])

    pickle.dump(carriers,open('/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/200k/carriers_90rec15fp10_10.pkl','wb'))
    # filtered_carriers = ef.filter_results(tops,carriers,2,3)
    glm = ef.SingleGLM()
    glm.set_cov(cov.pandas,context)
    pvs = ef.test_pipe(carriers,df,fam_ext,glm)
    pickle.dump(pvs,open('pvs_90rec15fp10_10_height.pkl','wb'))

    cov = ef.Covariates()
    cov.add_from_file('/sc/arion/projects/ipm2/for_roohy/UKBB.age_sex.pcs.sorted.covar')
    # cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/hdl_chol.ukbb.mapped.pheno')
    cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/chol.ukbb.mapped.pheno')

    glm = ef.SingleGLM()
    glm.set_cov(cov.pandas,context)
    pvs = ef.test_pipe(carriers,df,fam_ext,glm)
    pickle.dump(pvs,open('pvs_90rec15fp10_10_ldl.pkl','wb'))