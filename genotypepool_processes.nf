#! /usr/bin/env nextflow
nextflow.enable.dsl=2

process generatepoolvcf {
    container "${projectDir}/singularity/bcftools.sif"
    publishDir "${params.outdir}/poolvcf/", mode: 'copy', pattern: "*pool.vcf.gz*"

    input:
        val(chromosome)
        path(samplefile)

    output:
        tuple val(chromosome), path("${chromosome}.pool.vcf.gz"), path("${chromosome}.pool.vcf.gz.csi")
    
    script:
    """
    poolvcf="${chromosome}.pool.vcf.gz"

    bcftools view --threads ${params.threads} --types snps --force-samples -S ${samplefile} -Oz -o ${chromosome}.pool.vcf.gz \
    \$(compgen -f ${params.refvcfdir}/*${chromosome}*.vcf.gz | grep ".gz\$")
    bcftools index -c ${chromosome}.pool.vcf.gz
    """
}

process extractwc {
    container "${projectDir}/singularity/strandseq_Rtools.sif"
    publishDir "${params.outdir}", mode: 'copy', pattern: "wc_regions.txt"

    input:
        path(bamdir)

    output:
        path("wc_regions.txt")

    script:
    """
    #!/usr/bin/env Rscript
    wcfile <- paste0("${params.outdir}", "/wc_regions.txt")
    chromslist <- strsplit("${params.chromosomes}", ",")[[1]]

    if ( ! file.exists(wcfile)) {
        suppressPackageStartupMessages(library(breakpointR))
        suppressPackageStartupMessages(library(GenomicRanges))

        breakpointr(inputfolder="${bamdir}", outputfolder="${params.outdir}/BPR_output/", pairedEndReads=as.logical("${params.paired}"), 
            numCPU=${params.threads}, windowsize=175, binMethod="reads", background=0.15, 
            chromosomes=chromslist, maskRegions="${projectDir}/blacklist.highdepth.centromeres.bed")

        exportRegions(datapath=paste0("${params.outdir}", "/BPR_output/data"), file="wc_regions.txt", 
                    collapseInversions=FALSE, minRegionSize=10000, state="wc")
    } else {
        file.copy(wcfile, "./")
    }
    """
}

process mergeBam {
    container "${projectDir}/singularity/bowtie2_samtools_bedtools.sif"
    publishDir "${params.outdir}/mergedbam/", mode: 'copy', pattern: "*merge.bam*"

    input:
        val(chromosome)
        path(bamdir)
    
    output:
        tuple val(chromosome), path("${chromosome}.strandseq.merge.bam"), path("${chromosome}.strandseq.merge.bam.bai"), emit: mergedbam

    script:
    """
    mergedbamdir=${params.outdir}/mergedbam/
    strandmergebam="${chromosome}.strandseq.merge.bam"

    if [ -f \$mergedbamdir/\$strandmergebam ] && [ -f \$mergedbamdir/\$strandmergebam".bai" ]; then
        ln -s \$mergedbamdir/\$strandmergebam* ./

    else
        samtools merge -@ ${params.threads} -R ${chromosome} \$strandmergebam ${bamdir}/*bam 
        samtools index \$strandmergebam
    fi
    """
}

process snpsMergeBam {
    container "${projectDir}/singularity/bcftools.sif"
    publishDir "${params.outdir}/mergedsnps/", mode: 'copy', pattern: "*.vcf.gz*"
    
    input:
        tuple val(chromosome), path(bam), path(bamindex)
    
    output:
        tuple val(chromosome), path("${chromosome}.strandseq.merge.vcf.gz"), path("${chromosome}.strandseq.merge.vcf.gz.csi"), emit: mergevcf

    script:
    """
    mergedsnpsdir=${params.outdir}/mergedsnps/
    strandmergevcf="${chromosome}.strandseq.merge.vcf.gz"

    if [ -f \$mergedsnpsdir/\$strandmergevcf ] && [ -f \$mergedsnpsdir/\$strandmergevcf".csi" ]; then
        ln -s \$mergedsnpsdir/\$strandmergevcf* ./

    else
        bcftools mpileup --threads ${params.threads} -f ${params.ref} -Ou ${bam} | \
        bcftools call --threads ${params.threads} --skip-variants indels -mv -Ou | \
        bcftools view --types snps --genotype het --include 'INFO/DP>=10' -Oz -o \$strandmergevcf
        bcftools index -c \$strandmergevcf
    fi
    """
}

process genotypePool {
    container "${projectDir}/singularity/strandseq_Rtools.sif"
    publishDir "${params.outdir}/genotyperesults/", mode: 'copy', pattern: "*.tsv*"
    
    input:
        tuple val(chromosome), path(vcf), path(poolvcf)
        path(wcregion)
        path(bamdir)
    
    output:
        path("poolgenotype_${chromosome}.tsv"), emit: results

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(StrandPhaseR))
    source("${projectDir}/scripts/genotypeStrandScells.R")

    chrom=as.character("${chromosome}")
    chrnamedlist=list("${poolvcf}")
    names(chrnamedlist) <- chrom

    genotypes_df <- genotypeStrandScells(inputfolder="${bamdir}",
                                        strandS.vcf="${vcf}",
                                        popul.vcf.list=chrnamedlist,
                                        wc.regions="${wcregion}",
                                        chromosomes=chrom,
                                        min.snv.cov=2,
                                        max.snv.per.chr=5000)

    write.table(genotypes_df, paste0('poolgenotype_', chrom, '.tsv'), sep="\t", quote=FALSE)
    """

}

process assignSingleCell {
    container "${projectDir}/singularity/py310_viz.sif"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*CellLineAssign*"
    input:
        val(gtList)

    output:
        path('*CellLineAssign.tsv')

    script:
    """
    python ${projectDir}/scripts/cellLineStats.py '${gtList}'
    """
}