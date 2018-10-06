## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for GLASS-WG pipeline
## Authors: Floris Barthel, Nauman Javed
## A majority of the code used in this workflow was copied from Nauman Javed's
## GitHub repo: https://github.com/naumanjaved/fingerprint_maps
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

import subprocess
import os
import itertools as it
import sys
import argparse
import traceback
import time
import vcf

## Chromosomes
CHR = [str(j+1) for j in range(22)]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge haplotype map files from individual chromosomes
## See `README_phase3_callset_20150220` for more information on this dataset
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule haplotype_map_merge:
    input:
        expand("results/haplotype-map/final/{chr}.filtered.map", chr=CHR)
    output:
        "data/ref/fingerprint.filtered.map"
    conda:
        "../envs/haplotype.yaml"
    message:
        "Concatnating fingerprint maps"
    shell:
        "cat {config[haplotype_map][header]} {input} > {output}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This rule downloads the 1000GP variant callsets in VCF
## See `README_phase3_callset_20150220` for more information on this dataset
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule download_vcf:
    output:
        "data/1000GP-VCF/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        "data/1000GP-VCF/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
    params:
        dir = "data/1000GP-VCF"
    conda:
        "../envs/haplotype.yaml"
    message:
        "Downloading 1000 Genomes Phase 3 VCFs(hg19)\n"
        "Chromosome: {wildcards.chr}"
    shell:
        "cd {params.dir} && \
            wget -r -np -nH --cut-dirs 4 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz && \
            wget -r -np -nH --cut-dirs 4 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This rule downloads the 1000GP haplotypes
## See `https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html` for more information
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule download_map:
    output:
        "data/1000GP-map/genetic_map_chr{chr}_combined_b37.txt"
    params:
        dir = "data/1000GP-map"
    conda:
        "../envs/haplotype.yaml"
    message:
        "Downloading 1000 Genomes Phase 3 Recombination maps(hg19)\n"
        "Chromosome: {wildcards.chr}"
    shell:
        "cd {params.dir} && \
            wget -r -np -nH --cut-dirs 3 http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/genetic_map_chr{wildcards.chr}_combined_b37.txt"

rule extract_similar_SNPs:
    '''
    Writes to a new file the SNPs in an input VCF with population specific
    minor allele fractions(MAFs) that do not differ by more than SIM. For
    example,this function would reject a SNP with EUR_AF - AMR_AF  = 0.15 if
    SIM = 0.1.
    '''
    input:
        vcf = "data/1000GP-VCF/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        tbi = "data/1000GP-VCF/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
    output:
        "results/haplotype-map/snps/chr_{chr}-common-SNPs.list"
    message:
        "Creating list of SNPs with similar MAFs across populations..."
    run:
        SIM = config["haplotype_map"]["similarity"]
        vcf_reader = vcf.Reader(filename=input["vcf"])  # Loop over VCF records
        SNPs = []
        for record in vcf_reader:
            info = record.INFO
            if record.INFO['VT'][0] != 'SNP':  # Extract only SNPs
                continue

            SNP =str(record.ID)

            if SNP == '.' or not isinstance(SNP, str):
                continue  # Skip over SNPs with no ID

            if SNP in SNPs:  # Skip duplicates
                continue

            values = []  # List to hold population allele frequencies(AFs)
            differences = []  # List to hold pairwise AF differences

            for field in info.keys():  # Load population AFs
                if '_AF' in field:
                    values.append(info[field][0])

            differences = [abs(y-x) for x, y in it.combinations(values, 2)]
            if max(differences) > SIM or len(values) != 5 or SNP == ".":
                continue  # If similarity threshold is exceeded, skip SNP

            SNPs.append(SNP)

        ## Print SNPs to file
        with open(str(output), 'a') as similar_SNPs: 
            for snp in SNPs:
                similar_SNPs.write(snp + '\n')
                

rule create_VCFs:
    '''
    Extract SNPs from each chromosome that are:
    1) biallelic
    2) have population specific MAF differences no greater than 10% as
       determined by being listed in chr_i-common-SNPs.list
    3) are phased
    '''
    input:
        snplist = "results/haplotype-map/snps/chr_{chr}-common-SNPs.list",
        vcf = "data/1000GP-VCF/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    output:
        "results/haplotype-map/vcf/chr_{chr}.recode.vcf"
    params:
    	prefix = "results/haplotype-map/vcf/chr_{chr}"
    conda:
        "../envs/haplotype.yaml"
    message:
        "Creating filtered VCFs."
    log:
        "logs/haplotype-map/create_VCFs.chr{chr}.log"
    shell:
    	"vcftools --gzvcf {input.vcf} \
    		--maf {config[haplotype_map][min_maf]} \
    		--phased \
    		--remove-indels \
    		--min-alleles 2 \
    		--max-alleles 2 \
    		--recode \
    		--recode-INFO-all \
    		--snps {input.snplist} \
    		--out {params.prefix} \
    		> {log} 2>&1"
    
rule sort_VCF:
    '''
    Sorts VCF and removes duplicate IDs
    '''
    input:
        "results/haplotype-map/vcf/chr_{chr}.recode.vcf"
    output:
        vcf = "results/haplotype-map/vcf/{chr}.recode_u.vcf",
        dup = "results/haplotype-map/vcf/{chr}.duplicate"
    message:
        "Sorting VCFs..."
    run:
        VCF_file = str(input)
        
        # Sort VCF and record duplicates
        command = "awk -F'\t' '!($1 ~ /#/)' " + VCF_file \
                  + " | awk -F'\t' '{print $3}' | sort | uniq -d > " \
                  + output["dup"]
        subprocess.call(command, shell=True)

        duplicates = open(output["dup"], 'r')
        snplist = []  # Load duplicate IDs into snplist

        with duplicates as d:
            for snp in d:
                snplist.append(snp.rstrip('\n'))

        vcf_file = open(VCF_file, 'r')
        
        # Rewrite unique SNPs to new_VCF_file
        new_VCF_file = output
        newvcf = open(new_VCF_file, 'w')

        for line in vcf_file.readlines():

            if '#' in line.split('\t')[0]:  # Write header
                newvcf.write(line)
            else:
                if line.split('\t')[2] not in snplist:
                    newvcf.write(line)  # Write unique SNPs to new VCF

        newvcf.close()

rule create_PLINK_binary:
    '''
    Creates PLINK binary files from input VCF for use with
    LDSC script(https://github.com/bulik/ldsc).
    '''
    input:
        vcf = "results/haplotype-map/vcf/{chr}.recode_u.vcf",
        map = expand("data/1000GP-map/genetic_map_chr{chr}_combined_b37.txt", chr=CHR)
    output:
        "results/haplotype-map/plink/{chr}.bed"
    log:
        "logs/haplotype-map/create_PLINK_binary.chr{chr}.log"
    params:
        recomb_file = lambda wildcards: "data/1000GP-map/genetic_map_chr{chr}_combined_b37.txt* {chr}".format(chr = wildcards.chr)
    message:    
        "Creating PLINK binary files..."
    conda:
        "../envs/haplotype.yaml"
    shell:
        "plink \
            --debug \
            --vcf {input.vcf} \
            --cm-map {params.recomb_file} \
            --make-bed \
            --out results/haplotype-map/tmp/{wildcards.chr} \
            > {log} 2>&1"

rule LD_score:
    '''
    Calculates LDscore of all variants for each chromosome
    using LDSC script(https://github.com/bulik/ldsc)
    '''
    input:
        "results/haplotype-map/plink/{chr}.bed"
    output:
        "results/haplotype-map/ldscore/LD-{chr}.l2.ldscore.gz"
    params:
        bfile = "results/haplotype-map/plink/{chr}",
        script = lambda wildcards: config["haplotype_map"]["ld_script"],
        window = lambda wildcards: config["haplotype_map"]["ld_score_window_autosome"],
        prefix = "results/haplotype-map/ldscore/LD-{chr}"
    message:
        "Calculating LDScores..."
    log:
        "logs/haplotype-map/LD_score.chr{chr}.log"
    conda:
        "../envs/haplotype.yaml"
    shell:
        "python {params.script} \
            --bfile {params.bfile} \
            --ld-wind-cm {params.window} \
            --out {params.prefix} \
            --l2 --yes-really \
            > {log} 2>&1"

rule order:
    '''
    Here we reassign LDscore_i_new = (1.0 - LD_score_i / max(LD_score)).
    This gives the highest scoring SNP the lowest p value, and thus the
    greatest priority when clumping.
    '''
    input:
        "results/haplotype-map/ldscore/LD-{chr}.l2.ldscore.gz"
    output:
        "results/haplotype-map/order/{chr}.p"
    params:
        basename = "results/haplotype-map/ldscore/LD-{chr}.l2.ldscore"
    message:    
        "Creating PLINK association files..."
    run:
        association_file_name = params["basename"]
        
        # unzip LD score files
        subprocess.call("gzip -d " + input, shell=True)
        assoc = open(association_file_name, 'r')

        # write a new association file to rank SNPs by LDscore
        p_file = open(output, 'w')

        p_file.write("SNP" + "\t" + 'P\n')
        p_dictionary = {}  # Create a SNP:LDscore dictionary
        size_dict = {}

        for line in assoc.readlines()[1:]:
            SNP = line.split('\t')[1]
            LD = float(line.split('\t')[3].rstrip('\n'))
            p_dictionary[SNP] = LD
            assoc.close()
        
        # Calculate the max LDscore for the chromosome
        max_LD = float(max(p_dictionary.values()))
        
        '''
        For each SNP in the sorted dictionary, map the LDscore
        to LD_score_new_i = (1.0 - LD_score_i / max(LD_score))
        Add 1e-6 to set the minimum "signifigance value"
        '''
        for SNP in sorted(p_dictionary, key=p_dictionary.get, reverse=True):
            p_val = (1.0 - (float(p_dictionary[SNP])/(max_LD))) + 0.000001
            if p_val < 1.0:
                p_file.write(SNP + '\t' + str(p_val) + '\n')
        del p_dictionary
        p_file.close() 

rule prune:
    '''
    Take in all SNPs found on the chromosome and prune them to a set of
    independent variants. PLINK looks at all pairs of variants within an X
    variant window. If a pair has an R^2 correlation greater than some
    threshold, then one SNP within this pair is removed. After pruning within
    this window, PLINK slides along the chromosome by some number of SNPs.
    '''
    input:
        "results/haplotype-map/plink/{chr}.bed"
    output:
        "results/haplotype-map/prune/{chr}.prune.in"
    params:
        bfile = "results/haplotype-map/plink/{chr}",
        prefix = "results/haplotype-map/prune/{chr}"
    message:        
        "Pruning SNPs..."
    conda:
        "../envs/haplotype.yaml"
    log:
        "logs/haplotype-map/prune.chr{chr}.log"
    shell:
        "if [ \"{wildcards.chr}\" == \"X\" ]; \
        then \
            plink --bfile {params.bfile} \
                --indep-pairwise \
                {config[haplotype_map][prune_window]} \
                {config[haplotype_map][prune_slide]} \
                {config[haplotype_map][prune_cutoff]} \
                --r \
                --ld-xchr 1 \
                --out {params.prefix} \
                > {log} 2>&1; \
        else \
            plink --bfile {params.bfile} \
                --indep-pairwise \
                {config[haplotype_map][prune_window]} \
                {config[haplotype_map][prune_slide]} \
                {config[haplotype_map][prune_cutoff]} \
                --r \
                --out {params.prefix} \
                > {log} 2>&1; \
        fi"
    
rule LD_separate:
    '''
    Read in the association files created by order function
    and separate it out into independent(pruned.in) SNPs and
    dependent(pruned.out) SNPs.
    '''
    input:
        ordered = "results/haplotype-map/order/{chr}.p",
        pruned = "results/haplotype-map/prune/{chr}.prune.in"
    output:
        independent = "results/haplotype-map/separate/{chr}_in.p",
        dependent = "results/haplotype-map/separate/{chr}_out.p"
    message:    
        "Separating out LDScores into dependent and independent SNP files..."
    run:
        LD_file = open(input["ordered"], 'r')
        prune_file = open(input["pruned"], 'r')
        
        in_file = open(output["independent"], 'w')
        in_file.write('SNP' + '\t' + 'P\n')
        
        out_file = open(output["dependent"], 'w')
        out_file.write('SNP' + '\t' + 'P\n')
        
        out_SNPs = []
        index_SNPs = []

        with prune_file as prune:
            for k, line in enumerate(prune):
                index_SNPs.append(line.rstrip('\n'))

        with LD_file as LD:
            for k, line in enumerate(LD):
                if k > 0:
                    SNP = line.rstrip('\n').split('\t')[0]
                    if SNP in index_SNPs:
                        in_file.write(line)
                    else:
                        out_file.write(line)

        in_file.close()
        out_file.close()
    
rule clump:
    '''
    Clumps dependent variants(_out.p) around independent variants(_in.p) in a
    greedy manner. Each SNP "clumped" around an index variant should be within
    a certain window size(specified with clump-kb parameter) and have an r^2
    correlation of atleast "cutoff"(specified with --clump-r2) with the index
    variant. p1 and p2 represent the upper thresholds of signifigance, above
    which SNPs should be excluded. Here they are set to 1 to include all SNPs.
    '''
    input:
        independent = "results/haplotype-map/separate/{chr}_in.p",
        dependent = "results/haplotype-map/separate/{chr}_out.p"
    output:
        "results/haplotype-map/clump/{chr}.clumped",
        "results/haplotype-map/clump/{chr}.ld"
    params:
        bfile = "results/haplotype-map/plink/{chr}",
        outdir = "results/haplotype-map/clump/{chr}"
    message:    
        "Clumping SNPs..."
    conda:
        "../envs/haplotype.yaml"
    log:
        "logs/haplotype-map/clump.chr{chr}.log"
    shell:
        "plink --bfile {params.bfile} \
            --clump {input.independent} {input.dependent} \
            --clump-index-first \
            --clump-p1 1.0 \
            --clump-p2 1.0 \
            --clump-r2 {config[haplotype_map][clump_cutoff]} \
            --clump-kb {config[haplotype_map][max_distance_clump]} \
            --out {params.outdir} \
            > {log} 2>&1"
    
rule reformat_clumps:
    '''
    Script to read in clumps and write them into a haplotype map file readable
    by picard. This file format has recently been deprecated in favor of a VCF
    format file so this should be updated at some point.
    '''
    input:
        clumped = "results/haplotype-map/clump/{chr}.clumped",
        vcf = "results/haplotype-map/vcf/{chr}.recode_u.vcf"
    output:
        map = "results/haplotype-map/reformat/{chr}.map",
        sorted = "results/haplotype-map/reformat/{chr}_sorted.clumped"
    message:
        "Building map file..."
    run:
        # Sort the clumps by base pair position of index variant
        sorting_command = "tail -n+2 " + input["clumped"] \
                          + " | sort -k4,4 > " + output["sorted"]
        subprocess.check_call(sorting_command, shell=True)

        block_dict = {}  # Block dictionary with variant:index variant pairs
        anchor_file = open(output["sorted"], 'r')
        anchors = []
        for line in anchor_file.readlines()[2:]:
            anchor = line.rstrip('\n').split()[2]
            anchors.append(anchor)
            block_dict[anchor] = ""
            snps = line.rstrip('\n').split()[11].split(',')
            for snp in snps:
                if snp.rstrip('(w)') not in anchors:
                    variant = snp.split('(')[0]
                    block_dict[variant] = anchor
        anchor_file.close()

        old_map = vcf.Reader(filename=input["vcf"])
        new_map = open(output["map"], 'w')
        
        '''
        parse each line in old vcf file to extract the SNP name, base
        pair position, minor and major alleles, MAF and look up its index
        variant in the block dictionary defined above. If the variant IS
        an index variant then write its index variant as ''
        '''
        
        for record in old_map:
            if record.ID in block_dict.keys():
                newline = str(record.CHROM) + '\t' \
                          + str(record.POS) + '\t' \
                          + record.ID + '\t' \
                          + str(record.REF[0]) + '\t' \
                          + str(record.ALT[0]) + '\t' \
                          + str(record.INFO['AF'][0]) + '\t' \
                          + block_dict[record.ID] \
                          + '\t' + "" + '\t' + '\n'
                new_map.write(newline)
        new_map.close()
    
rule detect_negative_LD:
    '''
    Reads in data from clumping procedure and writes out pairs of SNPs in the
    newly created map file that are in high negative LD that give a high r^2
    correlation.
    '''
    input:
        map = "results/haplotype-map/reformat/{chr}.map",
        ld = "results/haplotype-map/clump/{chr}.ld"
    output:
        "results/haplotype-map/negld/{chr}.negLD"
    message:    
        "Detecting negative LD..."
    run:
        r_file = open(input["ld"], 'r')
        r_dict = {}  # Create a dictionary of SNPs:r^2 pairs
        
        with r_file as r:
            for k, line in enumerate(r):
                if k > 0:
                    r_dict[(line.split()[2], line.split()[5])] = \
                            line.split()[6].rstrip('\n')
        
        '''
        Loop over the newly created map file and for each
        pair of SNP/anchor SNP, determine whether the
        the pair's r^2 correlation is negative by checking
        its value in the dictionary above. If it is, write
        the pair to a new ".negLD" file
        '''

        map_file = open(input["map"], 'r')
        neg_LD_file_name = output
        with map_file as maps:
            for num, line in enumerate(maps):
                if len(line.split('\t')) > 6:
                    snp1 = line.split('\t')[2]
                    snp2 = line.split('\t')[6]
                    if (snp1, snp2) in r_dict.keys():
                        if float(r_dict[(snp1, snp2)]) < 0.0:
                            with open(neg_LD_file_name, 'a') as negLDfile:
                                negLDfile.write(snp1 + '\t' + snp2 + '\t'
                                                + r_dict[(snp1, snp2)] + '\n')
                                del r_dict[(snp1, snp2)]
                    if (snp2, snp1) in r_dict.keys():
                        if float(r_dict[(snp2, snp1)]) < 0.0:
                            with open(neg_LD_file_name, 'a') as negLDfile:
                                negLDfile.write(snp1 + '\t' + snp2 + '\t'
                                                + r_dict[(snp2, snp1)] + '\n')
                                del r_dict[(snp2, snp1)]
    
rule switch_alleles:
    '''
    For each chromosomal map file, read in the pairs of SNPs found to be in
    negative LD and switch the major and minor alleles in the map file.
    Keep original map file and write the new map file to *.filtered.map
    '''
    input:
        map = "results/haplotype-map/reformat/{chr}.map",
        negld = "results/haplotype-map/negld/{chr}.negLD"
    output:
        "results/haplotype-map/final/{chr}.filtered.map"
    message:    
        "Switching alleles..."
    run:
        # If no pairs were in negative linkage, copy to output folder and skip
        if not os.path.isfile(input["negld"]):
            copy_command = "cp " + input["map"] + " " + output
            subprocess.check_call(copy_command, shell=True)
            return None

        tab = '\t'
        negfile = open(input["negld"], 'r')
        neglist = []
        
        with negfile as neg:
            for line in neg:
                neglist.append(line.split()[0])
        
        newmapfile =open(output, 'w')
        oldmapfile = open(input["map"], 'r')
        
        with oldmapfile as old:
            for k, line in enumerate(old):
                chrom =line.split('\t')[0]
                bp = line.split('\t')[1]
                maf = line.split('\t')[5]
                snp = line.split('\t')[2]
                minor = line.split('\t')[4]
                maj = line.split('\t')[3]
                anch = line.split('\t')[6]
                # Switch alleles if found in the list of SNPs in negative LD
                if snp in neglist:
                    newline = [chrom, bp, snp, minor, maj, maf, anch, "", "\n"]
                    newmapfile.write(tab.join(newline))
                else:
                    newmapfile.write(line)
        newmapfile.close()

## END ##