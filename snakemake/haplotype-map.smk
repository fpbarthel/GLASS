
import map_building_functions as build

CHR = [str(j+1) for j in range(22)]

rule all:
    input: "output/fingerprint.filtered.map"

rule merge:
    input:
        header = "header",
        chrs = expand("output/{chr}.filtered.map", chr=CHR)
    output:
        "output/fingerprint.filtered.map"
    message:
        "Concatnating fingerprint maps"
    shell:
        "cat {input.header} {input.chrs} > {output}"


INTERMEDIATE_DIR = "intermediates/"
LD_SCRIPT = "/projects/barthf/GLASS-WG/sandbox/ldsc/ldsc.py"
RECOMB_DIR = "/projects/barthf/GLASS-WG/sandbox/1000GP_Phase3/"
FULL_PATH = "/projects/barthf/GLASS-WG/sandbox/fingerprint_maps/"
SIMILARITY = 0.10
MIN_MAF = 0.10
LD_SCORE_WINDOW_AUTOSOME = 1.0
PRUNE_WINDOW = 1500
PRUNE_SLIDE = 5
PRUNE_CUTOFF = 0.1
CLUMP_CUTOFF = 0.9
MAX_DISTANCE_CLUMP = 1500
  
rule extract_similar_SNPs:
    input:
        "/projects/barthf/GLASS-WG/sandbox/VCFs/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    output:
        "intermediates/chr_{chr}-common-SNPs.list"
    message:
        "Creating list of SNPs with similar MAFs across populations..."
    run:
        build.extract_similar_SNPs(wildcards["chr"], input, INTERMEDIATE_DIR, SIMILARITY)

rule create_VCFs:
    input:
        snplist = "intermediates/chr_{chr}-common-SNPs.list",
        vcf = "/projects/barthf/GLASS-WG/sandbox/VCFs/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    output:
        "intermediates/chr_{chr}.recode.vcf"
    message:
        "Creating filtered VCFs."
    run:
        build.create_VCFs(wildcards["chr"], input["vcf"], INTERMEDIATE_DIR, MIN_MAF)
    
rule sort_VCF:
    input:
        "intermediates/chr_{chr}.recode.vcf"
    output:
        "intermediates/{chr}.recode_u.vcf"
    message:
        "Sorting VCFs..."
    run:
        build.sort_VCF(wildcards["chr"], INTERMEDIATE_DIR)

rule create_PLINK_binary:
    input:
        "intermediates/{chr}.recode_u.vcf"
    output:
        "intermediates/{chr}.bed"
    log:
        "logs/create_PLINK_binary.chr{chr}.log"
    params:
        recomb_file = lambda wildcards: "{dir}genetic_map_chr{chr}_combined_b37.txt* {chr}".format(dir=RECOMB_DIR, chr = wildcards.chr)
    message:    
        "Creating PLINK binary files..."
    shell:
        "/projects/barthf/opt/miniconda3/envs/ldsc/bin/plink \
            --debug \
            --vcf {input} \
            --cm-map {params.recomb_file} \
            --make-bed \
            --out intermediates/{wildcards.chr} \
            > {log} 2>&1"

rule LD_score:
    input:
        "intermediates/{chr}.bed"
    output:
        "intermediates/LD-{chr}.l2.ldscore.gz"
    params:
        bfile = "intermediates/{chr}",
        script = lambda wildcards: LD_SCRIPT,
        window = lambda wildcards: LD_SCORE_WINDOW_AUTOSOME,
        prefix = "intermediates/LD-{chr}"
    message:
        "Calculating LDScores..."
    log:
        "logs/LD_score.chr{chr}.log"
    shell:
        "source activate ldsc; "
        "python {params.script} \
            --bfile {params.bfile} \
            --ld-wind-cm {params.window} \
            --out {params.prefix} \
            --l2 --yes-really \
            > {log} 2>&1"


rule order:
    input:
        "intermediates/LD-{chr}.l2.ldscore.gz"
    output:
        "intermediates/{chr}.p"
    message:    
        "Creating PLINK association files..."
    run:
        build.order(wildcards["chr"], INTERMEDIATE_DIR)

rule prune:
    input:
        "intermediates/{chr}.bed"
    output:
        "intermediates/{chr}.prune.in"
    message:        
        "Pruning SNPs..."
    run:
        build.prune(wildcards["chr"], INTERMEDIATE_DIR, PRUNE_WINDOW, PRUNE_SLIDE, PRUNE_CUTOFF)
    
rule LD_separate:
    input:
        "intermediates/{chr}.p",
        "intermediates/{chr}.prune.in"
    output:
        "intermediates/{chr}_in.p",
        "intermediates/{chr}_out.p"
    message:    
        "Separating out LDScores into dependent and independent SNP files..."
    run:
        build.LD_separate(wildcards["chr"], INTERMEDIATE_DIR)
    
rule clump:
    input:
        "intermediates/{chr}_in.p",
        "intermediates/{chr}_out.p"
    output:
        "intermediates/{chr}.clumped"
    message:    
        "Clumping SNPs..."
    run:
        build.clump(wildcards["chr"], INTERMEDIATE_DIR, CLUMP_CUTOFF, MAX_DISTANCE_CLUMP)
    
rule reformat_clumps:
    input:
        "intermediates/{chr}.clumped"
    output:
        "intermediates/{chr}.map"
    message:
        "Building map file..."
    run:
        build.reformat_clumps(wildcards["chr"], INTERMEDIATE_DIR)
    
rule detect_negative_LD:
    input:
        "intermediates/{chr}.map"
    output:
        "intermediates/{chr}.negLD"
    message:    
        "Detecting negative LD..."
    run:
        build.detect_negative_LD(wildcards["chr"], INTERMEDIATE_DIR)
    
rule switch_alleles:
    input:
        map = "intermediates/{chr}.map",
        negld = "intermediates/{chr}.negLD"
    output:
        "output/{chr}.filtered.map"
    message:    
        "Switching alleles..."
    run:
        build.switch_alleles(wildcards["chr"], FULL_PATH, INTERMEDIATE_DIR)

## END ##