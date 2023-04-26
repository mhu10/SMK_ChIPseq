configfile: "config/config.yaml"

IDS, = glob_wildcards(config['INPUT_DIR']+'{id}.fastq.gz')                ###---  modify filename if needed
SAMPLES, = glob_wildcards(config['INPUT_DIR']+'{sample}_R1.fastq.gz')     ###---  modify filename if needed


rule all:
    input:
        expand("results/fastqc/{id}.html", id=IDS),
        expand("results/fastqc/{id}_fastqc.zip", id=IDS),
        expand("results/trimmomatic/{sample}_R1_trimmed.fastq", sample=SAMPLES),
        expand("results/trimmomatic/{sample}_R2_trimmed.fastq", sample=SAMPLES),
        expand("results/trimmomatic/{sample}_R1_unpaired.fastq", sample=SAMPLES),
        expand("results/trimmomatic/{sample}_R2_unpaired.fastq", sample=SAMPLES),
	    expand("results/Bowtie2/{sample}_trimmed.sam",sample=SAMPLES),
	    expand("results/Samtools/{sample}_trimmed.bam",sample=SAMPLES),
	    expand("results/Samtools/{sample}_trimmed_sort.bam",sample=SAMPLES),
	    expand("results/Samtools/{sample}_trimmed_sort_q10.bam",sample=SAMPLES),
	    expand("results/Picard/{sample}_trimmed_sort_q10_rmdup.bam",sample=SAMPLES),
	    expand("results/Picard/{sample}_trimmed_sort_q10_rmdup.txt",sample=SAMPLES),
	    expand("results/MACS2/{sample}/{sample}_peaks.xls",sample=SAMPLES),
	    expand("results/MACS2/{sample}/{sample}_peaks.narrowPeak",sample=SAMPLES),
	    expand("results/MACS2/{sample}/{sample}_summits.bed",sample=SAMPLES),
	    expand("results/Picard/{sample}_trimmed_sort_q10_rmdup.bam.bai",sample=SAMPLES),
        expand("results/MultiBamSummary/deeptools_readCounts.npz"),
        expand("results/MultiBamSummary/deeptools_readCounts.tab"),
	    expand("results/MultiBamSummary/deepTools_pcaplot.png"),
	    expand("results/MultiBamSummary/deeptools_pcaProfile.tab"),
        expand("results/MultiBamSummary/deepTools_heatmap_spearman.png"),
        expand("results/HOMER/merged/merged_peaks_CTL.narrowPeak"),
        expand("results/HOMER/merged/merged_peaks_TREATMENT.narrowPeak"),
        expand("results/HOMER/merged/merged_peaks_CTL.bed"),
        expand("results/HOMER/merged/merged_peaks_TREATMENT.bed")
    

rule fastqc:
    input:
        config['INPUT_DIR']+'{id}.fastq.gz',
    output:
        html="results/fastqc/{id}.html",
        zip="results/fastqc/{id}_fastqc.zip",


    params: "--quiet"
    log:
        "logs/fastqc/{id}.log"

    threads: 8
    wrapper:
        "v1.7.1/bio/fastqc"


rule trimmomatic_pe:
    input:
        r1=config['INPUT_DIR']+'{sample}_R1.fastq.gz',
        r2=config['INPUT_DIR']+'{sample}_R2.fastq.gz'
    output:
        r1="results/trimmomatic/{sample}_R1_trimmed.fastq",
        r2="results/trimmomatic/{sample}_R2_trimmed.fastq",
        r1_unpaired="results/trimmomatic/{sample}_R1_unpaired.fastq",
        r2_unpaired="results/trimmomatic/{sample}_R2_unpaired.fastq"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:TruSeq3-PE.fa:2:30:10","LEADING:3","TRAILING:3","SLIDINGWINDOW:4:15","MINLEN:36"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        8
    resources:
        mem_mb=1024
    wrapper:
        "v1.7.1/bio/trimmomatic/pe"


rule bowtie2_map:
    input:
        r1="results/trimmomatic/{sample}_R1_trimmed.fastq",
        r2="results/trimmomatic/{sample}_R2_trimmed.fastq",
    params:
        bowtieindex=config['BOWTIE2INDEX']
    output:
        "results/Bowtie2/{sample}_trimmed.sam"
    threads: 8
    conda:
        'envs/BowtieMapping.yaml'
    shell:
        "bowtie2 -p 8 --local --no-discordant --no-mixed "
        "-x {params.bowtieindex}"
        "-1 {input.r1} "
        "-2 {input.r2} "
        "-S {output} "

rule samtools_bam:
    input:
        "results/Bowtie2/{sample}_trimmed.sam"
    output:
        "results/Samtools/{sample}_trimmed.bam"
    threads: 8
    conda:
        'envs/SamtoolsBAM.yaml'
    shell:
        "samtools view -bS {input} > {output} "


rule samtools_sort:
    input:
        "results/Samtools/{sample}_trimmed.bam"
    output:
        "results/Samtools/{sample}_trimmed_sort.bam"
    threads: 8
    conda:
        'envs/SamtoolsBAM.yaml'
    shell:
        "samtools sort -@ 5 {input}  -o {output} "


rule samtools_filter:
    input:
        "results/Samtools/{sample}_trimmed_sort.bam"
    output:
        "results/Samtools/{sample}_trimmed_sort_q10.bam"
    threads: 8
    conda:
        'envs/SamtoolsBAM.yaml'
    shell:
        "samtools view -bh -q 10 {input} > {output} "

rule mark_duplicates:
    input:
        bams="results/Samtools/{sample}_trimmed_sort_q10.bam"
    output:
        bam="results/Picard/{sample}_trimmed_sort_q10_rmdup.bam",
        metrics="results/Picard/{sample}_trimmed_sort_q10_rmdup.txt"
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wrapper:
        "v1.20.0/bio/picard/markduplicates"



rule samtools_index:
    input:
        "results/Picard/{sample}_trimmed_sort_q10_rmdup.bam",
    output:
        "results/Picard/{sample}_trimmed_sort_q10_rmdup.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 8  
    wrapper:
        "v1.25.0/bio/samtools/index"


rule MultiBamSummary:
    input:
        bams=expand("results/Picard/{sample}_trimmed_sort_q10_rmdup.bam",sample=SAMPLES),
        bais=expand("results/Picard/{sample}_trimmed_sort_q10_rmdup.bam.bai",sample=SAMPLES)
    output:
        countnpz = "results/MultiBamSummary/deeptools_readCounts.npz",
        counttab = "results/MultiBamSummary/deeptools_readCounts.tab"
    threads: 8
    conda:
        'envs/DeeptoolsQC.yaml'
    shell:
        "multiBamSummary bins "
        "--ignoreDuplicates -p 8 "
        "--bamfile {input.bams} "
        "-out {output.countnpz} "
        "--outRawCounts {output.counttab}"
        


rule MultiBamPCA:
    input:
        countnpz="results/MultiBamSummary/deeptools_readCounts.npz",
    output:
        "results/MultiBamSummary/deepTools_pcaplot.png",
        "results/MultiBamSummary/deeptools_pcaProfile.tab"
        
    threads: 8
    conda:
        'envs/DeeptoolsQC.yaml'
    shell:
        "plotPCA --corData {input.countnpz} "
        "--plotFile ./results/MultiBamSummary/deepTools_pcaplot.png "
        "-T 'PCA of read counts' "
        "--outFileNameData ./results/MultiBamSummary/deeptools_pcaProfile.tab "

rule MultiBamCorHeatmap:
    input:
        countnpz="results/MultiBamSummary/deeptools_readCounts.npz"
    output:
        "results/MultiBamSummary/deepTools_heatmap_spearman.png"
        
    threads: 8
    conda:
        'envs/DeeptoolsQC.yaml'
    shell:
        "plotCorrelation --corData {input.countnpz} "
        "--plotFile ./results/MultiBamSummary/deepTools_heatmap_spearman.png "
        "--corMethod spearman "
        "--whatToPlot heatmap "
        "--plotTitle 'Spearman Correlation of Read Counts' "
        "--plotNumbers "


rule MACS2_callpeak:
    input:
        "results/Picard/{sample}_trimmed_sort_q10_rmdup.bam"
    output:
        "results/MACS2/{sample}/{sample}_peaks.xls",
        "results/MACS2/{sample}/{sample}_peaks.narrowPeak",
        "results/MACS2/{sample}/{sample}_summits.bed",
    threads: 8
    conda:
        'envs/MACS2Callpeaks.yaml'
    shell:
        "macs2 callpeak -t {input} -n {wildcards.sample} --outdir results/MACS2/{wildcards.sample} -g hs -f BAM --keep-dup auto --bdg "


rule homer_mergePeaks_ctl:
    input:
        # input peak files
        expand("results/MACS2/{sample}/{sample}_peaks.narrowPeak",sample=config['HOMER_MERGE_CTL'])
    output:
        "results/HOMER/merged/merged_peaks_CTL.narrowPeak"
    params:
        extra="-d given"  # optional params, see homer manual
    log:
        "logs/mergePeaks/merged_peaks.log"
    wrapper:
        "v1.25.0/bio/homer/mergePeaks"


rule homer_mergePeaks_treatment:
    input:
        # input peak files
        expand("results/MACS2/{sample}/{sample}_peaks.narrowPeak",sample=config['HOMER_MERGE_TREATMENT'])
    output:
        "results/HOMER/merged/merged_peaks_TREATMENT.narrowPeak"
    params:
        extra="-d given"  # optional params, see homer manual
    log:
        "logs/mergePeaks/merged_peaks.log"
    wrapper:
        "v1.25.0/bio/homer/mergePeaks"


rule r_PeaksToBed:
    input:
        peaks_ctl="results/HOMER/merged/merged_peaks_CTL.narrowPeak",
        peaks_treatment="results/HOMER/merged/merged_peaks_TREATMENT.narrowPeak"
    output:
        bed_ctl="results/HOMER/merged/merged_peaks_CTL.bed",
        bed_treatment="results/HOMER/merged/merged_peaks_TREATMENT.bed"
    conda:
        "envs/r_PeaksToBed.yaml"
    script:
        "scripts/ConvertPeaks.R"
