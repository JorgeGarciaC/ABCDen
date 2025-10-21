#!/usr/bin/env python
# coding: utf-8

import sys
import allel
import re
import numpy as np

def Calculate_ZnS_python(vcf_file, init_pos=0, final_pos=210000):
	region_vcf='4'+':'+str(init_pos)+'-'+str(final_pos)
	vcf = allel.read_vcf(vcf_file, 
               ['variants/CHROM', 'variants/POS', 'calldata/GT', 'samples'],
			   region=region_vcf,tabix = None)
	Samples = vcf['samples']
	gt_py = allel.GenotypeArray(vcf['calldata/GT'])
	POS = vcf['variants/POS']

# Caculating ZnS
	def subset_gt(gt_py,population):
	    idx = list(map(lambda x: bool(re.match(population, x)), Samples))
	    gt_pop = gt_py.subset(sel1=idx)
	    is_seg = gt_pop.count_alleles().is_segregating() 
	    gt_s = gt_pop.compress(is_seg, axis=0)
	    return gt_s, POS[is_seg]

	gt_CHB, POS_CHB = subset_gt(gt_py,"CHB")    
	Seg_CHB = gt_CHB.shape[0]
	r2_CHB = allel.rogers_huff_r(gt_CHB.to_n_alt(fill=-1))**2
	ZnS_CHB = (np.nansum(r2_CHB)*2)/(Seg_CHB*(Seg_CHB-1))

	gt_EUR, POS_EUR = subset_gt(gt_py,"EUR")
	Seg_EUR = gt_EUR.shape[0]
	r2_EUR = allel.rogers_huff_r(gt_EUR.to_n_alt(fill=-1))**2
	ZnS_EUR = (np.nansum(r2_EUR)*2)/(Seg_EUR*(Seg_EUR-1))

	gt_Eafr, POS_Eafr = subset_gt(gt_py,"Eafr")
	Seg_Eafr = gt_Eafr.shape[0]
	r2_EAFR = allel.rogers_huff_r(gt_Eafr.to_n_alt(fill=-1))**2
	ZnS_EAFR = (np.nansum(r2_EAFR)*2)/(Seg_Eafr*(Seg_Eafr-1))
	
	gt_DENI, POS_DENI = subset_gt(gt_py,"DENI")
	S_DENI = gt_DENI.shape[0]
	r2_DENI = allel.rogers_huff_r(gt_DENI.to_n_alt(fill=-1))**2
	ZnS_DENI = (np.nansum(r2_DENI)*2)/(S_DENI*(S_DENI-1))
	
	NEA = np.array(list(map(lambda x: bool(re.match("NEA",x)),Samples)))
	VIN = np.array(list(map(lambda x: bool(re.match("VIN",x)),Samples)))
	gt_raw = gt_py.subset(sel1=NEA+VIN) # index of neanderthal samples
	is_seg = gt_raw.count_alleles().is_segregating() 
	gt_NEAN = gt_raw.compress(is_seg, axis=0)
	S_NEAN = gt_NEAN.shape[0]
	r2_NEAN = allel.rogers_huff_r(gt_NEAN.to_n_alt(fill=-1))**2
	ZnS_NEAN = (np.nansum(r2_NEAN)*2)/(S_NEAN*(S_NEAN-1))
	
	return(ZnS_CHB,ZnS_EUR,ZnS_EAFR, ZnS_DENI, ZnS_NEAN)
