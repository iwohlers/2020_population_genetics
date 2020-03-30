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

rule bcftools_stats:
    input: "{path}/{filename}.vcf.gz",
           "{path}/{filename}.vcf.gz.tbi"
    output: "{path}/{filename}.bcftools_stats"
    conda: "envs/bcftools.yaml"
    shell: "bcftools stats {input[0]} > {output}"

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
           "-Djava.io.tmpdir=\"/scratch/tmp\" " + \
           "-XX:ParallelGCThreads=4 " + \
           "-Xmx250g " + \
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
    input: "EGYPTGSA/controls_final.vcf.gz.tbi"


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
    input: "EGYPTGSAPSO/cases_final.vcf.gz.tbi"


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
    input: "EGYPTWGS/egyptians_final.vcf.gz.tbi"


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
    input: "HOLAZARIDIS2016/HumanOriginsPublic2068_final.vcf.gz.tbi"


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
    input: "BUSBY2020/BusbyWorldwidePopulations_final.vcf.gz.tbi"


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
    input: "1000G/1000G_final.vcf.gz.tbi"


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
    input: "FERNANDES2019/AP_IRAN_clean_final.vcf.gz.tbi"


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
    shell: "vcftools --gzvcf {input} " + \
                    "--remove-filtered-all " + \
                    "--recode " + \
                    "--stdout " + \
                    ">{output} 2>{log}"

rule preprocess_scott:
    input: "SCOTT2016/ciliopathies_exomes_2569_final.vcf.gz.tbi"


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
                   "hgdp_wgs.20190516.full.chr{wildcards.x}.vcf.gz " +\
                   ">/dev/null 2>&1"

# Downloading the meta data
rule download_bergstroem_metadata:
    output: "BERGSTROEM2020/hgdp_wgs.20190516.metadata.txt"
    shell: "wget -P BERGSTROEM2020 ftp://ngs.sanger.ac.uk/production/hgdp/" +\
                   "hgdp_wgs.20190516/metadata/" +\
                   "hgdp_wgs.20190516.metadata.txt"    

# Downloading all relevant files
rule download_bergstroem_all:
    input: expand("BERGSTROEM2020/hgdp_wgs.20190516.full.chr{x}.vcf.gz", \
                  x=[str(num) for num in range(1,23)]+["X","Y"]), \
           "BERGSTROEM2020/README.data-access.hgdp_wgs.20190516.txt", \
           "BERGSTROEM2020/hgdp_wgs.20190516.metadata.txt"

# Concatenate the vcf file from several chromosomes
rule concatenate_chr_vcfs_bergstroem:
    input: expand("BERGSTROEM2020/hgdp_wgs.20190516.full.chr{x}.vcf.gz", \
                  x=[str(num) for num in range(1,23)])
    output: "BERGSTROEM2020/hgdp_wgs.20190516.full_hg38.vcf.gz"
    conda: "envs/vcftools.yaml"
    shell: "vcf-concat {input} | bgzip > {output}"

rule preprocess_bergstroem:
    input: "BERGSTROEM2020/hgdp_wgs.20190516.full_final.vcf.gz.tbi"


################################################################################
############################ JOHN2018 data set #################################
################################################################################

rule download_john:
    output: "JOHN2018/Kuwaiti_Exomes_291.vcf.gz"
    shell: "wget -P JOHN2018 ftp://dgr.dasmaninstitute.org/Kuwaiti_Exomes_291.vcf.gz"

rule cp_unzip_replace_empty_john:
    input: "JOHN2018/Kuwaiti_Exomes_291.vcf.gz"
    output: "JOHN2018/Kuwaiti_Exomes_291_hg37.vcf"
    shell: "zcat {input} | sed 's/\t\t/\t.\t/g' > {output}"

rule preprocess_john:
    input: "JOHN2018/Kuwaiti_Exomes_291_final.vcf.gz.tbi"


################################################################################
######################## RODRIGUEZFLORES2016 data set ##########################
################################################################################

# For this dataset, we exclude variants with filter "LowQual", have to rename
# the chromosomes to include "chr" and perform variant lifting from 37 to 38
# Also, we may have to fix the VCF header

rule exclude_lowqual_add_header:
    input: "data/RODRIGUEZFLORES2016/QG108.Share/Vcf/QG108.Autosomal.gq30.dp10.vcf.gz",
           "data/RODRIGUEZFLORES2016/QG108.Share/Vcf/QG108.Autosomal.gq30.dp10.header.txt"
    output: "RODRIGUEZFLORES2016/QG108.Autosomal.gq30.dp10.vcf"
    shell: "cat {input[1]} > {output}; "
           "zcat {input[0]} | grep -v '##' | grep -v 'LowQual' >> {output} "

rule change_chrom_names_rodriguez:
    input: "RODRIGUEZFLORES2016/QG108.Autosomal.gq30.dp10.vcf",
           "liftover/change_chrom_names.txt"
    output: "RODRIGUEZFLORES2016/QG108.Autosomal.gq30.dp10_hg37.vcf"
    conda: "envs/bcftools.yaml"
    shell: "bcftools annotate --rename-chrs {input[1]} {input[0]} > {output}"

rule preprocess_rodriguez:
    input: "RODRIGUEZFLORES2016/QG108.Autosomal.gq30.dp10_final.vcf.gz.tbi"


################################################################################
############################ FAKHRO2018 data set ###############################
################################################################################

# For this exome sequencing dataset, we exclude variants with filter "LowQual", have to rename
# the chromosomes to include "chr" and perform variant lifting from 37 to 38
# Further, we exclude the individuals that have been wgs and are thus part of
# the RODRIGUEZFLORES2016 data set

rule exclude_lowqual:
    input: "data/FAKHRO2016/QG1005.Share/Integrated/Qatar.unrelated.all-1005.autosomal.integrated.SNPs.20160224.vcf.gz"
    output: "FAKHRO2016/Qatar.unrelated.all-1005.autosomal.integrated.SNPs.20160224.vcf"
    shell: "zcat {input} | grep '#' > {output}; "
           "zcat {input} | grep -v '#' | grep -v 'LowQual' >> {output} "

# We exclude all individuals that belong in fact to the RODRIGUESZFLORES206 data
# set, and all variants that after this individuals filtering have only 0/0 
# genotypes for all individuals (i.e. a minor allele count --mac of zero, i.e.
# --mac less than 1)
# This results in only 178477 variants.
# According to the Fakhro publication: "After exclusion of relatives and 
# application of batch-specific filters, an average of 4,045,064 SNPs were 
# observed per genome (n=88), an average of 15,382 SNPs were observed per 
# Exome51 Mb (n=853) and an average of 13,538 SNPs were observed per Exome38 Mb 
# (n=64). "
rule exclude_wgs_individuals:
    input: "FAKHRO2016/Qatar.unrelated.all-1005.autosomal.integrated.SNPs.20160224.vcf",
           "data/FAKHRO2016/QG1005.Share/Integrated/Qatar.unrelated.all-1005.autosomal.integrated.SNPs.20160224.exclude_individuals.txt"
    output: "FAKHRO2016/Qatar.unrelated.all-1005.autosomal.integrated.SNPs.20160224_exome.vcf"
    log: "FAKHRO2016/fakhro_exclude_wgs_individuals.log"
    conda: "envs/vcftools.yaml"
    shell: "vcftools --gzvcf {input[0]} " + \
                    "--remove {input[1]} " + \
                    "--mac 1 " + \
                    "--recode " + \
                    "--stdout " + \
                    ">{output} 2>{log}"

rule change_chrom_names_fakhro:
    input: "FAKHRO2016/Qatar.unrelated.all-1005.autosomal.integrated.SNPs.20160224_exome.vcf",
           "liftover/change_chrom_names.txt"
    output: "FAKHRO2016/Qatar.unrelated.all-1005.autosomal.integrated.SNPs.20160224_hg37.vcf"
    conda: "envs/bcftools.yaml"
    shell: "bcftools annotate --rename-chrs {input[1]} {input[0]} > {output}"

rule preprocess_fakhro:
    input: "FAKHRO2016/Qatar.unrelated.all-1005.autosomal.integrated.SNPs.20160224_final.vcf.gz.tbi"


################################################################################
#################### Compiling samples for analysis ############################
################################################################################

# Here, we compile the samples to be used for a specific admixture analyis 
# that is specified by a user-chosen analysis name
# The meta table columns by which to choose samples can be 
# SAMPLE, DATASET, WORLD_REGION, AFRICA_REGION, POPULATION, COUNTRY, 
# GENOTYPING_METHOD

rule compile_samples_for_analysis:
    input: "analysis_config/META_MASTER.txt",
           "analysis_config/config/{analysis}.txt"
    output: "admixture/{analysis}/sample_anno.txt"
    run:
        criteria = []
        with open(input[1],"r") as f_in:
            for line in f_in:
                if line[0] == '#' or line[:6] == 'SAMPLE':
                    continue
                criteria.append(line.strip("\n").split(","))
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                if line[:6] == "SAMPLE":
                    continue
                sample_meta = line.strip("\n").split("\t")
                for criterion in criteria:
                    ok = 0
                    for feature_num in range(7):
                        if criterion[feature_num] == '' or \
                           criterion[feature_num] == sample_meta[feature_num]:
                            ok += 1
                    if ok >= 7:
                        f_out.write(line)
                        break

rule dataset_list:
    input: "admixture/{analysis}/sample_anno.txt"
    output: "admixture/{analysis}/datasets.txt"
    shell: "cat {input} | cut -f 2 | sort | uniq > {output}"

DATASET_FILE = { \
    "1000G" : "1000G_final.vcf.gz", \
    "BERGSTROEM2020" : "hgdp_wgs.20190516.full_final.vcf.gz", \
    "BUSBY2020" : "BusbyWorldwidePopulations_final.vcf.gz", \
    "EGYPTGSA" : "controls_final.vcf.gz", \
    "EGYPTGSAPSO" : "cases_final.vcf.gz", \
    "EGYPTWGS" : "egyptians_final.vcf.gz", \
    "FAKHRO2016" : "Qatar.unrelated.all-1005.autosomal.integrated.SNPs.20160224_final.vcf.gz", \
    "FERNANDES2019" : "AP_IRAN_clean_final.vcf.gz", \
    "HOLAZARIDIS2016" : "HumanOriginsPublic2068_final.vcf.gz", \
    "JOHN2018" : "Kuwaiti_Exomes_291_final.vcf.gz", \
    "RODRIGUEZFLORES2016" : "QG108.Autosomal.gq30.dp10_final.vcf.gz", \
    "SCOTT2016" : "ciliopathies_exomes_2569_final.vcf.gz" \
}

rule symlinking_dataset_files:
     input: "admixture/{analysis}/datasets.txt"
     output: "admixture/{analysis}/data/symlinking.done"
     params: sym_prefix=lambda wildcards, output: output[0][:-15]
     run:
        shell("mkdir -p admixture/{wildcards.analysis}/data;")
        with open(input[0],"r") as f_in:
            for line in f_in:
                datasetname = line.strip("\n")
                target_file = DATASET_FILE[datasetname]
                # Symlink to variant files
                shell("ln -sf  ../../../"+datasetname+"/"+target_file+" "+ \
                      "admixture/"+wildcards.analysis+"/data/"+ \
                      datasetname+".vcf.gz")
                # Symlink to index files
                shell("ln -sf  ../../../"+datasetname+"/"+target_file+".tbi "+ \
                      "admixture/"+wildcards.analysis+"/data/"+ \
                      datasetname+".vcf.gz.tbi")
        shell("touch {output}")

# Intersecting files: bcftools isec options
# -c, --collapse snps|indels|both|all|some|none 
# -n, --nfiles [+-=]INT|~BITMAP
#    output positions present in this many (=), this many or more (+), this many 
#    or fewer (-), or the exact same (~) files 
rule intersecting_dataset_files:
    input:  "admixture/{analysis}/datasets.txt",
            "admixture/{analysis}/data/symlinking.done"
    output: "admixture/{analysis}/data/sites.txt",
            "admixture/{analysis}/data/0000.vcf.gz",
            "admixture/{analysis}/data/0001.vcf.gz",
    run: 
        # Intersecting all files
        # The bitmak used by bcftools isec should be all 1, because we want the
        # variants shared by all datasets; also, we don't want to merge anything
        # that is multiallelic
        bitmask = ""
        with open(input[0],"r") as f_in:
            for line in f_in:
                bitmask += "1" 
        # This generates files 0000.vcf,0001.vcf, etc
        shell("bcftools isec -c none " + \
                            "-n~"+bitmask+" " + \
                            "--output-type z " + \
                            "-p admixture/{wildcards.analysis}/data/ " + \
                            "admixture/{wildcards.analysis}/data/*.vcf.gz"
        )

# This rule gets only two data set files as input, but in fact merges up to 
# eleven (i.e. all available) data sets
# Merging file
#  -m, --merge snps|indels|both|all|none|id
#    The option controls what types of multiallelic records can be created: 
#         -m all    ..  SNP records can be merged with indel records
# In fact, if there are multialelic variants, then these should alrady have 
# identical ref and alt alleles after intersecting the files with option -c none
rule merging_datasets:
    input: "admixture/{analysis}/data/sites.txt",
           "admixture/{analysis}/data/0000.vcf.gz",
           "admixture/{analysis}/data/0001.vcf.gz"
    output: "admixture/{analysis}/{analysis}.vcf.gz"
    conda: "envs/bcftools.yaml"
    shell: "bcftools merge " + \
              "--merge all " + \
              "--output-type z " + \
              "--output {output} " + \
              "--threads 4 " + \
              "admixture/{wildcards.analysis}/data/00*.vcf.gz"

# Here, we apply actually several quite strict filtering criteria, in order to
# select variants that are not likely to have some technical artifacts, since
# we perform LD pruning anyway later, it is fine to exclude here all variants
# that may cause problems. We exclude:
# (i) indels (--remove-indels)
# (ii) variants with minor allele frequency less than parameter (--maf)
# (iii) variants with more than 5% missing genotypes (--max-missing)
# (iv) biallelic variants (--min-alleles 2 and --max-alleles 2)
# (v) significantly violating hardy weinberg disequilibrium (--hwe 0.000001)
# (vi) allow only autosomes, i.e. chromosomes 1-22
rule filter_for_admixture:
    input: "admixture/{analysis}/{analysis}.vcf.gz"
    output: "admixture/{analysis}/{analysis}_filtered_{maf}.vcf"
    log: "admixture/{analysis}/{analysis}_filtered_{maf}.log"
    wildcard_constraints: maf="([0-9]*[.])?[0-9]+"
    conda: "envs/vcftools.yaml"
    shell: "vcftools --gzvcf {input} " + \
                    "--remove-indels " + \
                    "--maf {wildcards.maf} " + \
                    "--max-missing 0.95 " + \
                    "--min-alleles 2 " + \
                    "--max-alleles 2 " + \
                    "--hwe 0.000001 " + \
                    "--chr chr1 --chr chr2 --chr chr3 --chr chr4 --chr chr5 "+ \
                    "--chr chr6 --chr chr7 --chr chr8 --chr chr9 "+ \
                    "--chr chr10 --chr chr11 --chr chr12 --chr chr13 " + \
                    "--chr chr14 --chr chr15 --chr chr16 --chr chr17 " + \
                    "--chr chr18 --chr chr19 --chr chr20 --chr chr21 " + \
                    "--chr chr22 "+ \
                    "--recode " + \
                    "--stdout " + \
                    "> {output} 2>{log}"

# Converting vcf files to plink binary format (bed/bim/fam) for admixture
# Sample ID conversion: 
# --double-id causes both family and individual IDs to be set to the sample ID
rule vcf_to_plink:
    input: "admixture/{analysis}/{analysis}_filtered_{maf}.vcf.gz"
    output: "admixture/{analysis}/{analysis}_filtered_{maf}.bed",
            "admixture/{analysis}/{analysis}_filtered_{maf}.bim",
            "admixture/{analysis}/{analysis}_filtered_{maf}.fam"
    log: "admixture/{analysis}/{analysis}_filtered_{maf}.vcf_to_plink.log"
    params: out_base=lambda wildcards, output: output[0][:-4]
    wildcard_constraints: maf="([0-9]*[.])?[0-9]+"
    conda: "envs/plink2.yaml"
    shell: "plink2 --vcf {input} " + \
                  "--double-id " + \
                  "--make-bed " + \
                  "--out {params.out_base} " + \
                  ">{log} 2>&1"

# LD prune the PLINK files; therefore, first make a list of SNPs in LD (and not 
# in LD)(i.e. to be removed or not)
# Parameters for indep-pairwise: [window size]<kb> [step size (variant ct)] 
# [VIF threshold]
# Explanation Plink website): the command above that specifies 50 5 0.5 would 
# a) consider a window of 50 SNPs, 
# b) calculate LD between each pair of SNPs in the window, 
# c) remove one of a pair of SNPs if the LD is greater than 0.5, 
# d) shift the window 5 SNPs forward and repeat the procedure
# Abraham 2014 used: 1000 10 0.02
# Anderson 2010 used: 50 5 0.2
# Wang 2009 used: 100 ? 0.2
# Fellay 2009 used: 1500 150 0.2 
# LD prune the PLINK files; therefore, first make a list of SNPs in LD (and not 
# in LD)(i.e. to be removed or not)
# Parameters for indep-pairwise: [window size]<kb> [step size (variant ct)] 
# [VIF threshold]
# Explanation Plink website): the command above that specifies 50 5 0.5 would 
# a) consider a window of 50 SNPs, 
# b) calculate LD between each pair of SNPs in the window, 
# c) remove one of a pair of SNPs if the LD is greater than 0.5, 
# d) shift the window 5 SNPs forward and repeat the procedure
# Abraham 2014 used: 1000 10 0.02
# Anderson 2010 used: 50 5 0.2
# Wang 2009 used: 100 ? 0.2
# Fellay 2009 used: 1500 150 0.2 
rule find_ld_pruned_snps:
    input: "admixture/{analysis}/{analysis}_filtered_{maf}.bed", 
           "admixture/{analysis}/{analysis}_filtered_{maf}.bim",
           "admixture/{analysis}/{analysis}_filtered_{maf}.fam"
    output: "admixture/{analysis}/{analysis}_filtered_{maf}.prune.in", 
            "admixture/{analysis}/{analysis}_filtered_{maf}.prune.out"
    log: "admixture/{analysis}/{analysis}_filtered_{maf}.find_ld_pruned_snps.log"
    params: in_base = lambda wildcards, input: input[0][:-4]
    conda: "envs/plink2.yaml"
    shell: "plink2 --bfile {params.in_base} " + \
                  "--indep-pairwise 1000 10 0.2 " + \
                  "--out {params.in_base} " + \
                  ">{log} 2>&1"

# Now exclude the pruned SNPs
rule exclude_ld_pruned_snps:
    input: "admixture/{analysis}/{analysis}_filtered_{maf}.bed", 
            "admixture/{analysis}/{analysis}_filtered_{maf}.bim",
            "admixture/{analysis}/{analysis}_filtered_{maf}.fam",
            "admixture/{analysis}/{analysis}_filtered_{maf}.prune.in"
    output: "admixture/{analysis}/{analysis}_filtered_{maf}_pruned.bed", 
            "admixture/{analysis}/{analysis}_filtered_{maf}_pruned.bim",
            "admixture/{analysis}/{analysis}_filtered_{maf}_pruned.fam"
    log: "admixture/{analysis}/{analysis}_filtered_{maf}.exclude_ld_pruned_snps.log"
    params: in_base = lambda wildcards, input: input[0][:-4],
            out_base = lambda wildcards, output: output[0][:-4]
    conda: "envs/plink2.yaml"
    shell: "plink2 --bfile {params.in_base} " + \
                  "--extract {input[3]} " + \
                  "--make-bed " + \
                  "--out {params.out_base} " + \
                  ">{log} 2>&1"


################################################################################
#################### Conducting admixture analysis #############################
################################################################################

# Getting the number of variants after intersecting the data sets,
# the number of variants after filtering and the number of variants after
# pruning, which is the number of variants going into admixture analysis.
rule num_variants_admixture:
    input: "admixture/{analysis}/{analysis}.vcf.gz",
           "admixture/{analysis}/{analysis}_filtered_{MAF}.vcf",
           "admixture/{analysis}/{analysis}_filtered_{MAF}_pruned.bim"
    output: "admixture/{analysis}/{analysis}_{MAF}_numvariants.txt"
    shell: "echo \"After intersection:\" > {output}; " + \
           "zcat {input[0]} | grep -v '#' | wc -l >> {output}; " + \
           "echo \"After filtering:\" >> {output}; " + \
           "cat {input[1]} | grep -v '#' | wc -l >> {output}; " + \ 
           "echo \"After pruning:\" >> {output}; " + \    
           "cat {input[2]} | wc -l >> {output}; "    

# Since admixture only produces the output files in the current directory, we
# go to the output dir and execute there
rule run_admixture:
    input: "admixture/{analysis}/{analysis}_filtered_{maf}_pruned.bed", 
    output: "admixture/{analysis}/{analysis}_filtered_{maf}_pruned.{K}.Q",
            "admixture/{analysis}/{analysis}_filtered_{maf}_pruned.{K}.P"
    log: "admixture/{analysis}/{analysis}_{maf}.{K}.log"
    wildcard_constraints: maf="([0-9]*[.])?[0-9]+", K="\d+"
    conda: "envs/admixture.yaml"
    shell: "cd admixture/{wildcards.analysis}; " + \
           "admixture  --seed=42 " + \
                       "-j48 " + \
                       "--cv=5 " + \
                       "../../{input[0]} {wildcards.K} > ../../{log}"

rule plot_admixture_q_values:
    input: "admixture/{analysis}/{analysis}_filtered_{maf}_pruned.{K}.Q",
           "admixture/{analysis}/{analysis}_filtered_{maf}_pruned.fam",
           "analysis_config/META_MASTER.txt"
    output: "admixture/{analysis}/{analysis}_{maf}.{K}.Q.pdf"
    conda: "envs/pophelper.yaml"
    wildcard_constraints: maf="([0-9]*[.])?[0-9]+", K="\d+"
    script: "scripts/plot_q_estimates.R"

rule compile_cv_values:
    input: expand("admixture/{{analysis}}/{{analysis}}_{{maf}}.{K}.log", \
                   K=range(1,16))
    output: "admixture/{analysis}/{analysis}_{maf}.cv"
    shell: "cat {input} | grep 'CV error' > {output}"

rule plot_admixture_pophelper:
    input: "analysis_config/META_MASTER.txt",
           "admixture/{analysis}/{analysis}_filtered_{maf}_pruned.fam",
           "admixture/{analysis}/{analysis}_{maf}.cv",
           "admixture/{analysis}/{analysis}_{maf}_numvariants.txt",
           expand("admixture/{{analysis}}/{{analysis}}_filtered_{{maf}}_pruned.{K}.Q", \
                   K=range(1,16))
    output: "admixture/{analysis}/{analysis}_{maf}.pophelper.pdf"
    params: out_path = lambda wildcards, output: "/".join(output[0].split("/")[:-1])+"/", \
            out_filename = lambda wildcards, output: output[0].split("/")[-1][:-4]
    conda: "envs/pophelper.yaml"
    wildcard_constraints: maf="([0-9]*[.])?[0-9]+", K="\d+"
    script: "scripts/pophelper.R"

rule admixture_all:
    input: expand("admixture/{analysis}/{analysis}_{maf}.pophelper.pdf", \
                   analysis = ["ADMIX_EGYPTGSA_EGYPTGSAPSO_BUSBY2020", \
                   "ADMIX_EGYPTGSA_EGYPTGSAPSO_1000G", \
                   "ADMIX_EGYPTGSA_EGYPTGSAPSO_BUSBY2020", \
                   "ADMIX_EGYPTGSA_EGYPTGSAPSO_EGYPTWGS", \
                   "ADMIX_EGYPTGSA_EGYPTGSAPSO_FAKHRO2016", \
                   "ADMIX_EGYPTGSA_EGYPTGSAPSO_FERNANDES2019", \
                   "ADMIX_EGYPTGSA_EGYPTGSAPSO_HOLAZARIDIS2016", \
                   "ADMIX_EGYPTGSA_EGYPTGSAPSO", \
                   "ADMIX_EGYPTGSA_EGYPTGSAPSO_SCOTT2016", \
                   "ADMIX_EGYPTWGS_1000G", \
                   "ADMIX_EGYPTWGS_BERGSTROEM2020"], \
                   maf=["0.00","0.05"])
    





