# We depart from merged vcfs of the ancient hominids
path_ancient = "/homes/users/jgarciac/projects/shared_data/OtherDatasets/Neanderthal/MERGE.ANCIENTS/"
path_modern = "/homes/users/jgarciac/projects/shared_data/1000G_2504_high_coverage/PANEL.38hg/GT/"


# Three positions
position = {"LEFT": "chr4:41896654-41966654", "CORE": "chr4:41966654-42036654", "RIGHT": "chr4:42036654-42106654"} 
names = ["LEFT","CORE","RIGHT"]

rule all:
    input:
        expand("obs.stats.chr{chrom}.{chunk}.tsv", chrom = ["4"], chunk = names)

# Important to fix problems when alt allele is not in a sample
rule fix_multi_allelic:
    input:
        anc = path_ancient + "4ne.chr{chrom}.nohmz.MASKED.lift.vcf.gz"
    output: 
        path_ancient + "4ne.chr{chrom}.nohmz.MASKED.norm.lift.vcf.gz"
    shell:
        "bcftools norm -m+any {input.anc} -Oz -o {output}; tabix {output}"

rule intersection:
    input: 
        anc = path_ancient + "4ne.chr{chrom}.nohmz.MASKED.norm.lift.vcf.gz",
        mod = path_modern + "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNVphased_panel.vcf.gz"
    output: 
        ancInter = "intersection.{chrom}.ancient.bed",
        modInter = "intersection.{chrom}.modern.bed"
    shell:
        "bedtools intersect -a {input.anc} -b {input.mod} -sorted > {output.ancInter};"
        "bedtools intersect -a {input.mod} -b {input.anc} -sorted > {output.modInter}"

rule fixALTallele:
    input:
        ancInter = "intersection.{chrom}.ancient.bed",
        modInter = "intersection.{chrom}.modern.bed"
    output:
        "new_ancientprevcf.{chrom}.txt"
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} FNR==NR {{pos[NR] = $5; next}} {{if ($5 == ".") $5 = pos[FNR]; print}}' {input.modInter} {input.ancInter} | sed 's/\//\|/g' > {output}
        """

rule genNewVcf:
    input: 
        preAncVcf="new_ancientprevcf.{chrom}.txt",
        preModVcf="intersection.{chrom}.modern.bed",
        anc = path_ancient + "4ne.chr{chrom}.nohmz.MASKED.norm.lift.vcf.gz",
        mod = path_modern + "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNVphased_panel.vcf.gz"
    output: 
        AncVcf = "anc.fix.{chrom}.vcf.gz",
        ModVcf = "mod.fix.{chrom}.vcf.gz"
    params: 
        aux= "fix.{chrom}.vcf"
    shell: 
        "bcftools view -h {input.anc} > header.{wildcards.chrom}.anc.txt;"
        "cat header.{wildcards.chrom}.anc.txt {input.preAncVcf} > anc.{params.aux};"
        "bgzip anc.{params.aux};"
        "tabix {output.AncVcf};"
        "bcftools view -h {input.mod} > header.{wildcards.chrom}.mod.txt;"
        "cat header.{wildcards.chrom}.mod.txt {input.preModVcf} > mod.{params.aux};"
        "bgzip mod.{params.aux};"
        "tabix {output.ModVcf}"

rule no_isec_merge:
    input:
        AncVcf = "anc.fix.{chrom}.vcf.gz",
        ModVcf = "mod.fix.{chrom}.vcf.gz"
    output:
        "merged.chr{chrom}.norm.vcf.gz"
    shell:
        "bcftools merge {input.AncVcf} {input.ModVcf} -Ou | bcftools view -m 2 -M 2 --type snps -Oz -o {output}; tabix {output}"
 
rule subset_nisc:
    input:
        "merged.chr{chrom}.norm.vcf.gz"
    output:
        "merged.chr{chrom}.norm.subsamples.vcf.gz"
    shell: 
        "bcftools view {input} -S subsample.txt --force-samples -Ou | bcftools view -m 2 -M 2 --type snps -Oz -o {output}; tabix {output}"

rule region:
    input:
        "merged.chr{chrom}.norm.subsamples.vcf.gz"
    output:
        "merged.chr{chrom}.{chunk}.norm.subset.vcf.gz"
    params:
        region=lambda wildcards: position[wildcards.chunk]  #"positions.get({chunk})" # large chr4:41851654-42151654
    shell:
        "bcftools view {input} -r {params.region} -Oz -o {output}; tabix {output}"

rule toHaps:
    input: 
        "merged.chr{chrom}.{chunk}.norm.subset.vcf.gz"
    output:
        haps=temp("merged.chr{chrom}.{chunk}.norm.subset.haps"),
        sample="merged.chr{chrom}.{chunk}.norm.subset.sample"
    params: 
        relate_path="~/scratch/relate1.1.x",
        vcf_name="merged.chr{chrom}.{chunk}.norm.subset"
    shell:
        "{params.relate_path}/bin/RelateFileFormats --mode ConvertFromVcf --haps {output.haps} --sample {output.sample} -i {params.vcf_name}"

rule flip:
    input: 
        haps="merged.chr{chrom}.{chunk}.norm.subset.haps",
        sample="merged.chr{chrom}.{chunk}.norm.subset.sample"
    output:
        "merged.chr{chrom}.{chunk}.norm.subset.flip.haps"
    params: 
        relate_path="~/scratch/relate1.1.x",
        ancestral_fa="~/projects/shared_data/1000G_2504_high_coverage/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_{chrom}.fa",
        out_name="merged.chr{chrom}.{chunk}.norm.subset.flip"
    shell:
        "{params.relate_path}/bin/RelateFileFormats --mode FlipHapsUsingAncestor "
        "--haps {input.haps} --sample {input.sample} "
        "--ancestor {params.ancestral_fa} -o {params.out_name}"        

rule region_only_modern:
    input:
        mod = path_modern + "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNVphased_panel.vcf.gz"
    output:
        "modern.chr{chrom}.{chunk}.subset.vcf.gz"
    params:
        region=lambda wildcards: position[wildcards.chunk] #"positions.get({chunk})"  #"chr4:41966654-42036654" # large chr4:41851654-42151654
    shell:
        "bcftools view {input.mod} -r {params.region} -S subsample.txt --force-samples -Ou | bcftools view - -m 2 -M 2 --type snps -Oz -o {output}; tabix {output}"

rule toHaps_mod:
    input: 
        "modern.chr{chrom}.{chunk}.subset.vcf.gz"
    output:
        haps="modern.chr{chrom}.{chunk}.subset.haps",
        sample="modern.chr{chrom}.{chunk}.subset.sample",
    params: 
        relate_path="~/scratch/relate1.1.x",
        vcf_name="modern.chr{chrom}.{chunk}.subset",
    shell:
        "{params.relate_path}/bin/RelateFileFormats --mode ConvertFromVcf --haps {output.haps} "
        "--sample {output.sample} -i {params.vcf_name}"

# Here I used Relate flippng utility as I find it really usefull
rule flipMod:
    input: 
        haps=temp("modern.chr{chrom}.{chunk}.subset.haps"),
        sample="modern.chr{chrom}.{chunk}.subset.sample"
    output: 
        flip_haps="modern.chr{chrom}.{chunk}.subset.flip.haps"
    params:
        relate_path="~/scratch/relate1.1.x",
        vcf_name="modern.chr{chrom}.{chunk}.subset",
        ancestral_fa="~/projects/shared_data/1000G_2504_high_coverage/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_{chrom}.fa",
    shell:
        "{params.relate_path}/bin/RelateFileFormats --mode FlipHapsUsingAncestor "
        "--haps {input.haps} --sample {input.sample} "
        "--ancestor {params.ancestral_fa} -o {params.vcf_name}.flip"

rule statsCalc:
    input: 
        flip_haps="modern.chr{chrom}.{chunk}.subset.flip.haps"
    output: 
        stats="obs.stats.chr{chrom}.{chunk}.tsv"
    params: 
        sing_img="/scratch/lab_ebosch/jgarciac/rstudio-singularity/rstudio-geo-R4.3.0_slendr-dev-parallel.sif"
    shell: 
        'singularity exec {params.sing_img} bash -c "Rscript Observed.stats.R 70000 {wildcards.chunk}"'	
