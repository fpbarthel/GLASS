## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for SomaticSeq Varscan2 and Mutect2 consensus calling
## Authors: Floris Barthel
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule somatiseq:
    input:
        vs2snp = "results/varscan2/fpfilter/{pair_barcode}.snp.Somatic.hc.final.vcf",
        vs2indel = "results/varscan2/vs2-filter/{pair_barcode}.indel.Somatic.hc.filter.vcf",
        mutect2 = "results/mutect2/final/{pair_barcode}.final.vcf",
        tumorbam = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        normalbam = lambda wildcards: "results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam".format(aliquot_barcode = manifest.getNormal(wildcards.pair_barcode))
    output:
        "results/somaticseq/{pair_barcode}/Consensus.sSNV.vcf"
        "results/somaticseq/{pair_barcode}/Consensus.sINDEL.vcf",
        "results/somaticseq/{pair_barcode}/Ensemble.sSNV.tsv",
        "results/somaticseq/{pair_barcode}/Ensemble.sINDEL.tsv"
    params:
    	outdir = "results/somaticseq/{pair_barcode}",
        mem = CLUSTER_META["somaticseq"]["mem"]
    threads:
        CLUSTER_META["somaticseq"]["ppn"]
    conda:
        "../envs/somaticseq.yaml"
    log:
        "logs/somaticseq/{pair_barcode}.log"
    benchmark:
        "benchmarks/somaticseq/{pair_barcode}.txt"
    message:
        "Running SomaticSeq consensus calling\n"
        "Pair: {wildcards.pair_barcode}"
    shell:
        "SomaticSeq.Wrapper.sh \
            --output-directory {params.outdir} \
            --genome-reference {config[reference_fasta]} \
            paired \
            --tumor-bam-file {input.tumorbam} \
			--normal-bam-file {input.normalbam} \
			--mutect2-vcf {input.mutect2} \
			--varscan-snv {input.vs2snp} \
			--varscan-indel {input.vs2indel} \
            > {log} 2>&1; "