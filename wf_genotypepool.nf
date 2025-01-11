#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETER CHECK
if ( ! params.samplefile ) {
    error("Error: --samplefile parameter not provided.")
}
if ( ! params.paired ) {
    error("Error: --paired (true/false) parameter not specified.")
}
if ( ! params.bamdir ) {
    error("Error: --bamdir parameter not provided.")
}
if ( ! params.ref ) {
    error("Error: --ref parameter not provided.")
}

// MESSAGE
if ( ! params.poolvcf ){
    println "Pool VCF not created."
}

// PRCOESSES
include {
    generatepoolvcf
    extractwc
    mergeBam
    snpsMergeBam
    genotypePool
    assignSingleCell
} from './genotypepool_processes.nf'

// HELP MESSAGE
if ( params.help ) {
    help = """wf_genotypepool.nf: To genotype individual single cell to 
             |
             |Required arguments:
             |  --samplefile  File with list of samples in pool.
             |  --paired      Paired end reads (true/false)
             |  --bamdir      Path to BAM (single-cell) directory. 
             |  --ref         Path to reference genome fasta. 
             |  --refvcfdir   Path to directory holding 1KG vcfs.
             | 
             |Optional arguments:
             |  --outdir      Output (Pool) Directory. [default: "./"]
             |  --poolvcf     Create pool VCF (Subset from larger VCF file). 
             |  --poolvcfdir  Absolute path to pool VCF. [default: ${params.poolvcfdir}]
             |  --chromosomes Comma-separated list of chromosomes to process.
             |                [default: ${params.chromosomes}]
             |  --threads     Number of threads. [default: ${params.threads}]
             |  --conda       Path to miniconda envs [default: ${params.conda}]""".stripMargin()

    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

// SET UP
log.info """\
    POOL GENOTYPE PIPELINE
    ===================================
    started at          : ${workflow.start}
    ---
    output directory    : ${params.outdir}
    bam directory       : ${params.bamdir}
    sample file         : ${params.samplefile}
    reference genome    : ${params.ref}
    paired-end reads    : ${params.paired}
    threads             : ${params.threads}
    """
    .stripIndent(true)



workflow {
    // Channels
    chrom_ch = Channel.fromList("${params.chromosomes}".split(",").toList())

    // Extract WC_regions
    wcregion_ch = extractwc("${params.bamdir}").collect()

    // Create EMBL Pool
    if ( ! params.poolvcf ){
        poolvcf_ch = Channel.fromPath("${params.poolvcfdir}/chr*.pool.vcf.gz").map{ file -> tuple(file.baseName.split("\\.")[0], file, "${file}.csi") }
    }
    else {
        samplefile_ch = Channel.fromPath("${params.samplefile}").collect()
        poolvcf_ch = generatepoolvcf(chrom_ch, samplefile_ch)
    }

    // Create merged BAM File
    mergebam_ch = mergeBam(chrom_ch)
    
    // Extract SNPs from merged BAM
    mergevcf_ch = snpsMergeBam(mergebam_ch)

    // Genotype Pool
    genotype_ch = mergevcf_ch.join(poolvcf_ch).map{ chrom, vcf, vcfidx, poolvcf, poolvcfidx -> tuple(chrom, vcf, poolvcf)}
    results_ch = genotypePool(genotype_ch, wcregion_ch)

    assignSingleCell(results_ch.collect())
} 