# kate:syntax python;
#
# This workflow performs population genetics analysis for the MENA region
# combining various data sets:
# 110 Egyptian individuals WGS
# ...

import gzip

rule all:
    input: ""


################################################################################
################# Some general helper files used more than once ################
################################################################################

# Changing chromosome names: without "chr" to with "chr"
rule file_for_change_chrom_names:
    output: "liftover/change_chrom_names.txt"
    shell: "echo \"1 chr1\" > {output}; " + \
           "echo \"2 chr2\" >> {output}; " + \
           "echo \"3 chr3\" >> {output}; " + \
           "echo \"4 chr4\" >> {output}; " + \
           "echo \"5 chr5\" >> {output}; " + \
           "echo \"6 chr6\" >> {output}; " + \
           "echo \"7 chr7\" >> {output}; " + \
           "echo \"8 chr8\" >> {output}; " + \
           "echo \"9 chr9\" >> {output}; " + \
           "echo \"10 chr10\" >> {output}; " + \
           "echo \"11 chr11\" >> {output}; " + \
           "echo \"12 chr12\" >> {output}; " + \
           "echo \"13 chr13\" >> {output}; " + \
           "echo \"14 chr14\" >> {output}; " + \
           "echo \"15 chr15\" >> {output}; " + \
           "echo \"16 chr16\" >> {output}; " + \
           "echo \"17 chr17\" >> {output}; " + \
           "echo \"18 chr18\" >> {output}; " + \
           "echo \"19 chr19\" >> {output}; " + \
           "echo \"20 chr20\" >> {output}; " + \
           "echo \"21 chr21\" >> {output}; " + \
           "echo \"22 chr22\" >> {output}"

rule compress_for_bcftools:
    input: "{DATASET}/{vcf_file}.vcf"
    output: "{DATASET}/{vcf_file}.vcf.gz"
    conda: "envs/bcftools.yaml"
    shell: "cat {input} | bgzip > {output}"

rule index_for_bcftools:
    input: "{DATASET}/{vcf_file}.vcf.gz"
    output: "{DATASET}/{vcf_file}.vcf.gz.tbi"
    conda: "envs/bcftools.yaml"
    shell: "tabix -p vcf {input}"

rule list_samplenames:
    input: "{DATASET}/{vcf_file}.vcf.gz",
           "{DATASET}/{vcf_file}.vcf.gz.tbi"
    output: "{DATASET}/{vcf_file}.samplenames"
    conda: "envs/bcftools.yaml"
    shell: "bcftools query --list-samples {input[0]} > {output}"

rule file_for_updating_samplenames:
    input: "{DATASET}/{vcf_file}.samplenames"
    output: "{DATASET}/{vcf_file}.updatesamplenames"
    run:
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                f_out.write(line.strip("\n")+"\t"+wildcards.DATASET+"_"+line)

rule update_samplenames:
    input: "{DATASET}/{vcf_file}_hg38.vcf.gz",
           "{DATASET}/{vcf_file}_hg38.vcf.gz.tbi",
           "{DATASET}/{vcf_file}_hg38.updatesamplenames"
    output: "{DATASET}/{vcf_file}_final.vcf.gz"
    conda: "envs/bcftools.yaml"
    shell: "bcftools reheader --samples {input[2]} {input[0]} > {output[0]}"    


################################################################################
############################## Variant liftover ################################
################################################################################

# Getting the chain file for liftover
rule get_liftover_chain_file_19To38:
    output: "liftover/hg19ToHg38.over.chain.gz"
    shell: "wget -P liftover ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

# Getting the chain file for liftover
rule get_liftover_chain_file_18To38:
    output: "liftover/hg18ToHg38.over.chain.gz"
    shell: "wget -P liftover ftp://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz"

# Getting the ucsc reference sequence for liftover
rule get_ucsc_reference_sequence:
    output: "liftover/hg38.fa.gz"
    shell: "wget -P liftover https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"

# Makeing a sequence dict needed by picard
rule make_dict:
    input: "liftover/hg38.fa.gz"
    output: "liftover/hg38.fa.gz.dict"
    conda: "envs/picard.yaml"
    shell: "picard CreateSequenceDictionary R={input} O={output}"    

# For liftover it is necessary to provide the java virtual machine with 
# additional memory, here to 280G: -Xmx280g
rule liftover_37To38:
    input: vcf="{dataset}/{filename}_hg37.vcf",
           chain="liftover/hg19ToHg38.over.chain.gz",
           dict="liftover/hg38.fa.gz.dict",
           ref="liftover/hg38.fa.gz"
    output: lifted="{dataset}/{filename}_hg38.vcf",
            rejected="{dataset}/{filename}_rejected.vcf",
    log: "{dataset}/{filename}.log"
    conda: "envs/picard.yaml"
    shell: "picard LiftoverVcf " + \
           "-Xmx280g " + \
           "I={input.vcf} " + \
           "O={output.lifted} " + \
           "CHAIN={input.chain} " + \
           "REJECT={output.rejected} " + \
           "RECOVER_SWAPPED_REF_ALT=true " + \
           "R={input.ref} > {log} 2>&1"

# For liftover it is necessary to provide the java virtual machine with 
# additional memory, here to 280G: -Xmx280g
rule liftover_36To38:
    input: vcf="{dataset}/{filename}_hg36.vcf",
           chain="liftover/hg18ToHg38.over.chain.gz",
           dict="liftover/hg38.fa.gz.dict",
           ref="liftover/hg38.fa.gz"
    output: lifted="{dataset}/{filename}_hg38.vcf",
            rejected="{dataset}/{filename}_rejected.vcf",
    log: "{dataset}/{filename}.log"
    conda: "envs/picard.yaml"
    shell: "picard LiftoverVcf " + \
           "-Xmx280g " + \
           "I={input.vcf} " + \
           "O={output.lifted} " + \
           "CHAIN={input.chain} " + \
           "REJECT={output.rejected} " + \
           "RECOVER_SWAPPED_REF_ALT=true " + \
           "R={input.ref} > {log} 2>&1"    


################################################################################
############################# EGYPTGSA data set ################################
################################################################################

rule change_chrom_names_egyptgsa:
    input: "data/unimputed/all.geno.mind.hetfail.king.pca.controls.vcf.gz",
           "liftover/change_chrom_names.txt"
    output: "EGYPTGSA/controls_hg37.vcf"
    conda: "envs/bcftools.yaml"
    shell: "bcftools annotate --rename-chrs {input[1]} {input[0]} > {output}"

rule preprocess_egyptgsa:
    input: "EGYPTGSA/controls_hg38.vcf.gz.tbi"


################################################################################
######################### EGYPTGSA psoriasis data set ##########################
################################################################################

rule change_chrom_names_egyptgsapso:
    input: "data/unimputed/all.geno.mind.hetfail.king.pca.cases.vcf.gz",
           "liftover/change_chrom_names.txt"
    output: "EGYPTGSAPSO/cases_hg37.vcf"
    conda: "envs/bcftools.yaml"
    shell: "bcftools annotate --rename-chrs {input[1]} {input[0]} > {output}"

rule preprocess_egyptgsapso:
    input: "EGYPTGSAPSO/cases_hg38.vcf.gz.tbi"


################################################################################
############################ EGYPTWGS data set ################################
################################################################################

rule egypt_wgs:
    input: "EGYPTWGS/egyptians_hg38.vcf"

# The variant calling by Matthias (data/vars.clean.vcf.gz) already contains
# the chromosome names with trailing "chr". Thus, here we just unzip
rule cp_and_unzip_egypt_wgs:
    input: "data/vars.clean.vcf.gz"
    output: "EGYPTWGS/egyptians_hg38.vcf"
    shell: "zcat {input} > {output}"

rule preprocess_egyptwgs:
    input: "EGYPTWGS/egyptians_hg38.vcf.gz.tbi"


################################################################################
####################### HOLAZARIDIS2016 data set ############################
################################################################################

# Getting the Near East genotype data from Lazaridis et al. Nature 2016
# (https://www.nature.com/articles/nature19310)
rule get_lazaridis_neareast_data:
    output: "HOLAZARIDIS2016/NearEastPublic.tar.gz"
    shell: "wget -P HOLAZARIDIS2016 https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/NearEastPublic.tar.gz"

rule extract_lazaridis_neareast_data:
    input: "HOLAZARIDIS2016/NearEastPublic.tar.gz"
    output: "HOLAZARIDIS2016/HumanOriginsPublic2068.geno",
            "HOLAZARIDIS2016/HumanOriginsPublic2068.snp",
            "HOLAZARIDIS2016/HumanOriginsPublic2068.ind",
            "HOLAZARIDIS2016/Permissions.txt"
    shell: "tar -xvzf {input} -C HOLAZARIDIS2016/"

rule write_convertf_par_file:
    input: "HOLAZARIDIS2016/HumanOriginsPublic2068.geno",
           "HOLAZARIDIS2016/HumanOriginsPublic2068.snp",
           "HOLAZARIDIS2016/HumanOriginsPublic2068.ind"
    output: "HOLAZARIDIS2016/par.PACKEDANCESTRYMAP.PED"
    shell: "echo 'genotypename:    {input[0]}' > {output}; " + \
           "echo 'snpname:         {input[1]}' >> {output}; " + \
           "echo 'indivname:       {input[2]}' >> {output}; " + \
           "echo 'outputformat:    PED' >> {output}; " + \
           "echo 'genotypeoutname: HOLAZARIDIS2016/HumanOriginsPublic2068.ped' >> {output}; " + \
           "echo 'snpoutname:      HOLAZARIDIS2016/HumanOriginsPublic2068.map' >> {output}; " + \
           "echo 'indivoutname:    HOLAZARIDIS2016/HumanOriginsPublic2068.pedind' >> {output}; "

# Convert from nativ ancestrymap to plink ped/map
rule convert_lazaridis:
    input: "HOLAZARIDIS2016/par.PACKEDANCESTRYMAP.PED",
           "HOLAZARIDIS2016/HumanOriginsPublic2068.geno",
           "HOLAZARIDIS2016/HumanOriginsPublic2068.ind",
           "HOLAZARIDIS2016/HumanOriginsPublic2068.snp"
    output: "HOLAZARIDIS2016/HumanOriginsPublic2068.ped",
            "HOLAZARIDIS2016/HumanOriginsPublic2068.map",
            "HOLAZARIDIS2016/HumanOriginsPublic2068.pedind"
    conda: "envs/eigensoft.yaml"
    shell: "convertf -p {input[0]}"

# Convert from plink ped/map to vcf
# vcf-iid: use within-family IDs fo sample IDs
rule lazaridis_ped2vcf:
    input: "HOLAZARIDIS2016/HumanOriginsPublic2068.ped",
           "HOLAZARIDIS2016/HumanOriginsPublic2068.map",
           "HOLAZARIDIS2016/HumanOriginsPublic2068.pedind"
    output: "HOLAZARIDIS2016/HumanOriginsPublic2068_wo_chr.vcf"
    log: "HOLAZARIDIS2016/lazaridis_ped2vcf.log"
    params: prefix=lambda wildcards, output: output[0][:-4]
    conda: "envs/plink2.yaml"
    shell: "plink2 --ped {input[0]} " + \
                  "--map {input[1]} " + \
                  "--recode vcf-iid " + \
                  "--real-ref-alleles " + \
                  "--alleleACGT " + \
                  "--out {params.prefix} > {log} 2>&1"

rule change_chrom_names_lazaridis:
    input: "HOLAZARIDIS2016/HumanOriginsPublic2068_wo_chr.vcf",
           "liftover/change_chrom_names.txt"
    output: "HOLAZARIDIS2016/HumanOriginsPublic2068_hg37.vcf"
    conda: "envs/bcftools.yaml"
    shell: "bcftools annotate --rename-chrs {input[1]} {input[0]} > {output}"

rule preprocess_lazaridis:
    input: "HOLAZARIDIS2016/HumanOriginsPublic2068_hg38.vcf.gz.tbi"


################################################################################
############################# BUSBY data set ###################################
################################################################################

# Getting Busby data
rule get_busby_data:
    output: "BUSBY2020/ckz9mtgrjj-3.zip"
    shell: "wget -P BUSBY2020/ https://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/ckz9mtgrjj-3.zip"

rule unzip_busby_data:
    input: "BUSBY2020/ckz9mtgrjj-3.zip"
    output: "BUSBY2020/BusbyWorldwide.corrected_ids.txt",
            "BUSBY2020/BusbyWorldwidePopulations_old.csv",
            "BUSBY2020/BusbyWorldwidePopulations.bed",
            "BUSBY2020/BusbyWorldwidePopulations.bim",
            "BUSBY2020/BusbyWorldwidePopulations.fam",
            "BUSBY2020/BusbyWorldwidePopulations_summary.corrected.xlsx"
    shell: "unzip {input} -d BUSBY2020"

# Some individuals (n=30) are duplicates within the Busby data set. These are 
# denoted by the author and they have the same individual ID, but different 
# family IDs. Those family IDs are in fact the population. Since family IDs
# are not used anymore in the VCF, we remove the duplicates as a first step.
# Also, after removing duplicates, a vcf file is produced
rule busby_bedbimfam2vcf:
    input: "BUSBY2020/BusbyWorldwidePopulations.bed",
           "BUSBY2020/BusbyWorldwidePopulations.bim",
           "BUSBY2020/BusbyWorldwidePopulations.fam",
           "data/busby_duplicates.txt"
    output: "BUSBY2020/BusbyWorldwidePopulations_wo_chr.vcf"
    log: "BUSBY2020/busby_bedbimfam2vcf.log"
    params: prefix_out=lambda wildcards, output: output[0][:-4]
    conda: "envs/plink2.yaml"
    shell: "plink2 --bed {input[0]} " + \
                  "--bim {input[1]} " + \
                  "--fam {input[2]} " + \
                  "--remove {input[3]} " + \
                  "--recode vcf-iid " + \
                  "--real-ref-alleles " + \
                  "--alleleACGT " + \
                  "--out {params.prefix_out} > {log} 2>&1"

rule change_chrom_names_busby:
    input: "BUSBY2020/BusbyWorldwidePopulations_wo_chr.vcf",
           "liftover/change_chrom_names.txt"
    output: "BUSBY2020/BusbyWorldwidePopulations_hg36.vcf"
    conda: "envs/bcftools.yaml"
    shell: "bcftools annotate --rename-chrs {input[1]} {input[0]} > {output}"

rule preprocess_busby:
    input: "BUSBY2020/BusbyWorldwidePopulations_hg38.vcf.gz.tbi"


################################################################################
############################# 1000G data set ###################################
################################################################################

# Downloading 1000 genomes data
rule download_1000g_genotypes:
    output: "1000G/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz"
    shell: "wget -P 1000G/ http://ftp.1000genomes.ebi.ac.uk/vol1/" + \
                                  "ftp/release/20130502/supporting/" + \
                                  "GRCh38_positions/" + \
                                  "ALL.chr{wildcards.chr}_GRCh38.genotypes.20170504.vcf.gz"

# Downloading 1000 genomes data (index)
rule download_1000g_genotypes_index:
    output: "1000G/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz.tbi"
    shell: "wget -P 1000G/ http://ftp.1000genomes.ebi.ac.uk/vol1/" + \
                                  "ftp/release/20130502/supporting/" + \
                                  "GRCh38_positions/" + \
                                  "ALL.chr{wildcards.chr}_GRCh38.genotypes.20170504.vcf.gz.tbi"

# Downloading 1000 genomes data (Readme)
rule download_1000g_genotypes_readme:
    output: "1000G/README_GRCh38_liftover_20170504.txt"
    shell: "wget -P 1000G/ http://ftp.1000genomes.ebi.ac.uk/vol1/" + \
                                  "ftp/release/20130502/supporting/" + \
                                  "GRCh38_positions/" + \
                                  "README_GRCh38_liftover_20170504.txt"

# Get the ped file which contains the population of the samples (and more info)
rule download_1000g_genotypes_ped:
    output: "1000G/integrated_call_samples_v2.20130502.ALL.ped"
    shell: "wget -P 1000G/ http://ftp.1000genomes.ebi.ac.uk/vol1/" + \
                                  "ftp/release/20130502/" + \
                                  "integrated_call_samples_v2.20130502.ALL.ped"

# This is downloading all chromosomes VCF files of 1000G, and two meta files
# Watch out: we already now use only chromosomes 1 to 22, thus no X,Y,MT,... 
rule download_1000g_genotypes_all:
    input: expand("1000G/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz", \
                   chr=[str(x) for x in range(1,23)]), \
           expand("1000G/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz.tbi", \
                   chr=[str(x) for x in range(1,23)]), \
           "1000G/README_GRCh38_liftover_20170504.txt", \
           "1000G/integrated_call_samples_v2.20130502.ALL.ped"

# Concatenate the vcf file from several chromosomes
rule concatenate_chr_vcfs_1000g:
    input: expand("1000G/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz", \
                   chr=[str(x) for x in range(1,23)])
    output: "1000G/1000G.vcf.gz"
    conda: "envs/vcftools.yaml"
    shell: "vcf-concat {input} | bgzip > {output}"

rule preprocess_1000g:
    input: "1000G/1000G.vcf.gz.tbi"


################################################################################
####################### FERNANDES2019 data set #################################
################################################################################

# Convert bed/bim/fam to vcf
rule fernandes_bedbimfam2vcf:
    input: "data/FERNANDES_2019_GRCh37/AP_IRAN_clean.bed",
           "data/FERNANDES_2019_GRCh37/AP_IRAN_clean.bim",
           "data/FERNANDES_2019_GRCh37/AP_IRAN_clean.fam"
    output: "FERNANDES2019/AP_IRAN_clean_wo_chr.vcf"
    params: prefix_out=lambda wildcards, output: output[0][:-4]
    log: "FERNANDES2019/fernandes_bedbimfam2vcf.log"
    conda: "envs/plink2.yaml"
    shell: "plink2 --bed {input[0]} " + \
                  "--bim {input[1]} " + \
                  "--fam {input[2]} " + \
                  "--recode vcf-iid " + \
                  "--real-ref-alleles " + \
                  "--alleleACGT " + \
                  "--out {params.prefix_out} > {log} 2>&1"

rule change_chrom_names_fernandes:
    input: "FERNANDES2019/AP_IRAN_clean_wo_chr.vcf",
           "liftover/change_chrom_names.txt"
    output: "FERNANDES2019/AP_IRAN_clean_hg37.vcf"
    conda: "envs/bcftools.yaml"
    shell: "bcftools annotate --rename-chrs {input[1]} {input[0]} > {output}"

rule preprocess_fernandes:
    input: "FERNANDES2019/AP_IRAN_clean_hg38.vcf.gz.tbi"


################################################################################
########################### SCOTT2016 data set #################################
################################################################################

# Filter to include only PASS variants
# --remove-filtered-all: Removes all sites with a FILTER flag other than PASS.
rule keep_pass_variants:
    input: "data/SCOTT_2016/74223/PhenoGenotypeFiles/"+ \
           "RootStudyConsentSet_phs000288.Ciliopathies_Exome.v2.p2.c1.GRU/"+\
           "GenotypeFiles/vcf/ciliopathies_exomes_2569.vcf.gz"
    output: "SCOTT2016/ciliopathies_exomes_2569_hg37.vcf"
    log: "SCOTT2016/scott_keeppass.log"
    conda: "envs/vcftools.yaml"
    params: prefix_out=lambda wildcards, output: output[0][:-4]
    shell: "vcftools --gzvcf {input} " + \
                    "--remove-filtered-all " + \
                    "--recode " + \
                    "--out {params.prefix_out} > {log} 2>&1"

rule preprocess_scott:
    input: "SCOTT2016/ciliopathies_exomes_2569_hg38.vcf.gz.tbi"


################################################################################
######################## BERGSTROEM2020 data set ###############################
################################################################################

# Downloading the README file
rule download_bergstroem_anno:
    output: "BERGSTROEM2020/README.data-access.hgdp_wgs.20190516.txt"
    shell: "wget -P BERGSTROEM2020 ftp://ngs.sanger.ac.uk/production/hgdp/" +\
                   "hgdp_wgs.20190516/README.data-access.hgdp_wgs.20190516.txt"

# Downloading the VCF files
rule download_bergstroem_vcf:
    output: "BERGSTROEM2020/hgdp_wgs.20190516.full.chr{x}.vcf.gz"
    shell: "wget -P BERGSTROEM2020 ftp://ngs.sanger.ac.uk/production/hgdp/" +\
                   "hgdp_wgs.20190516/" +\
                   "hgdp_wgs.20190516.full.chr{wildcards.x}.vcf.gz"

# Downloading the meta data
rule download_bergstroem_metadata:
    output: "BERGSTROEM2020/hgdp_wgs.20190516.metadata.txt"
    shell: "wget -P BERGSTROEM2020 ftp://ngs.sanger.ac.uk/production/hgdp/" +\
                   "hgdp_wgs.20190516/metadata/" +\
                   "hgdp_wgs.20190516.metadata.txt >/dev/null"    

# Downloading all relevant files
rule download_bergstroem_all:
    input: expand("BERGSTROEM2020/hgdp_wgs.20190516.full.chr{x}.vcf.gz", \
                  x=[str(num) for num in range(1,23)]+["X","Y"]), \
           "BERGSTROEM2020/README.data-access.hgdp_wgs.20190516.txt", \
           "BBERGSTROEM2020/hgdp_wgs.20190516.metadata.txt"

# Concatenate the vcf file from several chromosomes
rule concatenate_chr_vcfs_bergstroem:
    input: expand("BERGSTROEM2020/hgdp_wgs.20190516.full.chr{x}.vcf.gz", \
                  x=[str(num) for num in range(1,23)])
    output: "BERGSTROEM2020/hgdp_wgs.20190516.full_hg38.vcf.gz"
    conda: "envs/vcftools.yaml"
    shell: "vcf-concat {input} | bgzip > {output}"

rule preprocess_bergstroem:
    input: "BERGSTROEM2020/hgdp_wgs.20190516.full_hg38.vcf.gz.tbi"
