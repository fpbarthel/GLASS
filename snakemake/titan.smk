rule preparetitan:
    input:
        seg     = lambda wildcards: "results/cnv/mergedenoisedreadcounts/{aliquot_barcode}.denoisedCR.tsv".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode)),
        hets    = lambda wildcards: "results/cnv/modelsegments/{aliquot_barcode}/{aliquot_barcode}.hets.tsv".format(aliquot_barcode = manifest.getTumor(wildcards.pair_barcode))
    output:
        hets    = "results/cnv/titan/{pair_barcode}/{pair_barcode}.hets.tsv",
        seg     = "results/cnv/titan/{pair_barcode}/{pair_barcode}.seg"
    threads:
        CLUSTER_META["preparetitan"]["ppn"]
    log:
        "logs/cnv/preparetitan/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/preparetitan/{pair_barcode}.txt"
    message:
        "Prepare TitanCNA\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:"""
        cat {input.hets} | awk '/^[^@]/ {{ print $1,$2,$5,$3,$6,$4 }}' | tr ' ' '\\t' > {output.hets}
        cat {input.seg} | awk -F\\t 'BEGIN{{print"chr\\tstart\\tend\\tlog2_TNratio_corrected"}} /^[^@C]/ {{ print $0 }}' > {output.seg}
    """

rule titan:
    input:
        hets    = "results/cnv/titan/{pair_barcode}/{pair_barcode}.hets.tsv",
        seg     = "results/cnv/titan/{pair_barcode}/{pair_barcode}.seg"
    output:
        seg     = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.seg",
        segs    = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.segs.txt",
        titan   = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.titan.txt",
        params  = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.params.txt",
        cf      = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_CF.pdf",
        cna     = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_CNA.pdf",
        cnaseg  = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_CNASEG.pdf",
        loh     = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_LOH.pdf",
        lohseg  = "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}/{pair_barcode}_cluster0{cluster}_LOHSEG.pdf"
    params:
        outdir =        "results/cnv/titan/{pair_barcode}/ploidy{ploidy}/",
        alphaK =        lambda wildcards: 2500 if manifest.isExome(manifest.getTumor(wildcards.pair_barcode)) else 10000,
        alphaKHigh =    lambda wildcards: 2500 if manifest.isExome(manifest.getTumor(wildcards.pair_barcode)) else 10000,
        gender = 		lambda wildcards: manifest.getSex(manifest.getTumor(wildcards.pair_barcode)) if manifest.getSex(manifest.getTumor(wildcards.pair_barcode)) is not None else "NA"
    threads:
        CLUSTER_META["titan"]["ppn"]
    log:
        "logs/cnv/titan/{pair_barcode}_ploidy{ploidy}_cluster{cluster}.log"
    benchmark:
        "benchmarks/cnv/titan/{pair_barcode}_ploidy{ploidy}_cluster{cluster}.txt"
    message:
        "Run TitanCNA\n"
        "Pair ID: {wildcards.pair_barcode}\n"
        "Cluster: {wildcards.cluster}\n"
        "Ploidy: {wildcards.ploidy}"
    shell:"""
        set +u
        source activate R341
        set -u
        Rscript /projects/barthf/opt/TitanCNA/scripts/R_scripts/titanCNA.R \
            --genomeBuild hg19 \
            --id {wildcards.pair_barcode} \
            --hetFile {input.hets} \
            --cnFile {input.seg} \
            --numClusters {wildcards.cluster} \
            --numCores {threads} \
            --normal_0 0.5 \
            --ploidy_0 {wildcards.ploidy} \
            --alphaK {params.alphaK} \
            --alphaKHigh {params.alphaKHigh} \
            --minDepth 5 \
            --maxDepth 50000 \
            --gender {params.gender} \
            --estimatePloidy TRUE \
            --outDir {params.outdir} \
            --libdir /projects/barthf/opt/TitanCNA \
            > {log} 2>&1
    """

rule selecttitan:
    input:
        seg     = lambda wildcards: expand("results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.seg", ploidy = [2,3], cluster = [1,2,3], pair_barcode = wildcards.pair_barcode),
        segs    = lambda wildcards: expand("results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.segs.txt", ploidy = [2,3], cluster = [1,2,3], pair_barcode = wildcards.pair_barcode),
        titan   = lambda wildcards: expand("results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.titan.txt", ploidy = [2,3], cluster = [1,2,3], pair_barcode = wildcards.pair_barcode),
        params  = lambda wildcards: expand("results/cnv/titan/{pair_barcode}/ploidy{ploidy}/{pair_barcode}_cluster0{cluster}.params.txt", ploidy = [2,3], cluster = [1,2,3], pair_barcode = wildcards.pair_barcode)
    output:
        txt     = "results/cnv/titan/{pair_barcode}/{pair_barcode}.optimalClusters.txt"
    threads:
        CLUSTER_META["selecttitan"]["ppn"]
    log:
        "logs/cnv/selecttitan/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/selecttitan/{pair_barcode}.txt"
    message:
        "Select optimal TitanCNA cluster and ploidy\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:"""
        set +u
        source activate R341
        set -u
        Rscript /projects/barthf/opt/TitanCNA/scripts/R_scripts/selectSolution.R \
            --ploidyRun2=results/cnv/titan/{wildcards.pair_barcode}/ploidy2 \
            --ploidyRun3=results/cnv/titan/{wildcards.pair_barcode}/ploidy3 \
            --threshold=0.15 \
            --outFile {output.txt} \
            > {log} 2>&1
    """


rule finaltitan:
    input:
        txt     = "results/cnv/titan/{pair_barcode}/{pair_barcode}.optimalClusters.txt"
    output:
        segs    = "results/cnv/titanfinal/seg/{pair_barcode}.seg.txt",
        params  = "results/cnv/titanfinal/params/{pair_barcode}.params.txt",
        pdf     = "results/cnv/titanfinal/pdf/{pair_barcode}.pdf",
        igv		= "results/cnv/titanfinal/igv/{pair_barcode}.igv.seg"
    threads:
        CLUSTER_META["finaltitan"]["ppn"]
    log:
        "logs/cnv/finaltitan/{pair_barcode}.log"
    benchmark:
        "benchmarks/cnv/finaltitan/{pair_barcode}.txt"
    message:
        "Copy selected (final) TitanCNA results\n"
        "Pair ID: {wildcards.pair_barcode}"
    shell:"""
        TITANDIR=$(cat {input} | awk -F'\\t' '{{print $11}}' | sed -n 2p)
        cp "${{TITANDIR}}.segs.txt" {output.segs}
        cp "${{TITANDIR}}.seg" {output.igv}
        cp "${{TITANDIR}}.params.txt" {output.params}
        /opt/software/helix/ImageMagick/7.0.7-26/bin/montage $TITANDIR/*_CNA.pdf $TITANDIR/*_LOH.pdf $TITANDIR/*_CF.pdf -tile 1x3 -geometry +0+0 {output.pdf}
    """

## END ##