import sys
import allel
import re
import numpy as np
import pandas

def Calculate_ZnS_python(vcf_file,pop_labs,POP):
    vcf = allel.read_vcf(vcf_file, 
               ['variants/CHROM', 'variants/POS', 'calldata/GT', 'samples'])
    samples = vcf['samples']
    gt_py = allel.GenotypeArray(vcf['calldata/GT'])
    POS = vcf['variants/POS']
    panel = pandas.read_csv(pop_labs, sep=' ', names=['sample','pop'])
    panel = panel[panel['sample'].isin(samples)]
    samples_list = list(samples)
    samples_callset_index = [samples_list.index(s) for s in panel['sample']]
    panel['vcf_index'] = samples_callset_index
#    EUR_samples = panel[panel['pop'] == 'GBR'].index.values
#    Eafr_samples = panel[panel['pop'] == 'LWK'].index.values
#    CHB_samples = panel[panel['pop'] == 'CHB'].index.values
#    DEN_samples = panel[panel['pop'] == 'DEN'].index.values
#    NEA_samples = panel[panel['pop'] == 'NEA'].index.values
    samples = panel[panel['pop'] == POP].index.values
    
# Caculating ZnS
    def subset_gt(gt_py,idx):
        gt_pop = gt_py.subset(sel1=idx)
        is_seg = gt_pop.count_alleles().is_segregating() 
        gt_s = gt_pop.compress(is_seg, axis=0)
        return gt_s, POS[is_seg]

    gt_s, POS_s = subset_gt(gt_py,samples)    
    Seg_s = gt_s.shape[0]
    r2_s = allel.rogers_huff_r(gt_s.to_n_alt(fill=-1))**2
    ZnS_s = (np.nansum(r2_s)*2)/(Seg_s*(Seg_s-1))
    return(ZnS_s)
#    gt_CHB, POS_CHB = subset_gt(gt_py,CHB_samples)    
 #   Seg_CHB = gt_CHB.shape[0]
 #   r2_CHB = allel.rogers_huff_r(gt_CHB.to_n_alt(fill=-1))**2
 #   ZnS_CHB = (np.nansum(r2_CHB)*2)/(Seg_CHB*(Seg_CHB-1))

#    gt_EUR, POS_EUR = subset_gt(gt_py,EUR_samples)
  #  Seg_EUR = gt_EUR.shape[0]
   # r2_EUR = allel.rogers_huff_r(gt_EUR.to_n_alt(fill=-1))**2
    #ZnS_EUR = (np.nansum(r2_EUR)*2)/(Seg_EUR*(Seg_EUR-1))

    #gt_Eafr, POS_Eafr = subset_gt(gt_py,Eafr_samples)
   # Seg_Eafr = gt_Eafr.shape[0]
   # r2_EAFR = allel.rogers_huff_r(gt_Eafr.to_n_alt(fill=-1))**2
   # ZnS_EAFR = (np.nansum(r2_EAFR)*2)/(Seg_Eafr*(Seg_Eafr-1))
    
   # gt_DENI, POS_DENI = subset_gt(gt_py,DEN_samples)
   # S_DENI = gt_DENI.shape[0]
   # r2_DENI = allel.rogers_huff_r(gt_DENI.to_n_alt(fill=-1))**2
   # ZnS_DENI = (np.nansum(r2_DENI)*2)/(S_DENI*(S_DENI-1))
    
   # gt_NEAN, POS_NEAN = subset_gt(gt_py,NEA_samples)
   # S_NEAN = gt_NEAN.shape[0]
   # r2_NEAN = allel.rogers_huff_r(gt_NEAN.to_n_alt(fill=-1))**2
   # ZnS_NEAN = (np.nansum(r2_NEAN)*2)/(S_NEAN*(S_NEAN-1))
    
    return(ZnS_CHB,ZnS_EUR,ZnS_EAFR, ZnS_DENI, ZnS_NEAN)
