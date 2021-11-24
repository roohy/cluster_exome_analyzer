import pickle,sys
import extract_family as ef
import pyarrow.feather as feather
if __name__ == '__main__':
    # carriers = pickle.load(open('/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/200k/carriers_tops.pkl','rb'))
    carriers = pickle.load(open('/sc/arion/projects/ipm2/roohy/ukbb_ibd/exome/200k/carriers_90rec15fp10_10.pkl','rb'))
    cov = ef.Covariates()
    cov.add_from_file('/sc/arion/projects/ipm2/for_roohy/UKBB.age_sex.pcs.sorted.covar')
    # cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/hdl_chol.ukbb.mapped.pheno')
    # cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/standing_height.ukbb.mapped.pheno')
    cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/urae.ukbb_mapped.pheno')
    context = ef.ExomeContext('/sc/arion/projects/ipm2/roohy/ukbb_ibd/global_network/inon/output_chr1.id_list')
    fam_ext = ef.FamilyExtractor('/sc/arion/projects/ipm2/roohy/ukbb_ibd/umapping_25_i2_cls/')
    car_ext = ef.CarrierExtractor('/sc/arion/projects/ipm2/roohy/ukbb_ibd/global_network/inon/',context)
    tops = feather.read_feather( 'full_rec_no_fp_10_10.feather')
    # filtered_carriers = ef.filter_results(tops,carriers,2,3)
    glm = ef.SingleGLM()
    glm.set_cov(cov.pandas,context)
    pvs = ef.test_pipe(carriers,tops,fam_ext,glm)
    pickle.dump(pvs,open('pvs_90_15_10_10urae.pkl','wb'))


    cov = ef.Covariates()
    cov.add_from_file('/sc/arion/projects/ipm2/for_roohy/UKBB.age_sex.pcs.sorted.covar')
    # cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/hdl_chol.ukbb.mapped.pheno')
    # cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/chol.ukbb.mapped.pheno')
    cov.add_pheno('/sc/arion/projects/ipm2/for_roohy/creatinine.ukbb.mapped.pheno')
    glm = ef.SingleGLM()
    glm.set_cov(cov.pandas,context)
    pvs = ef.test_pipe(carriers,tops,fam_ext,glm)
    pickle.dump(pvs,open('pvs_90_15_10_10_creatinine.pkl','wb'))