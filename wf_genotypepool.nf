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
             |  --bamdir      Path to BAM directory single cell libraries. 
             |  --ref         Path to reference genome fasta. 
             |  --refvcfdir   Path to directory holding 1KG vcfs.
             | 
             |Optional arguments:
             |  --outdir      Output (Pool) Directory. [default: "./"]
             |  --poolvcfdir  Absolute path to pool VCF. [default: ${params.poolvcfdir}]
             |  --wcregion    File with WC regions of libraries with format: chr\tstart\tend\tlibrary. 
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
    Channel.fromList("${params.chromosomes}".split(",").toList()).set{ chroms }
    Channel.fromPath("${params.bamdir}").collect().set{ bamdir }

    // Extract WC_regions
    if ( params.wcregion ){
        Channel.fromPath("${params.wcregion}").set{ wcregion }
    }
    else {
        wcregion = extractwc(bamdir)
    }
    
    // Create EMBL Pool
    if ( params.poolvcfdir ){
        Channel.fromPath("${params.poolvcfdir}/chr*.pool.vcf.gz")
            .map{ file -> tuple(file.baseName.split("\\.")[0], file, "${file}.csi") }
            .set{ poolvcfs }
    }
    else {
        Channel.fromPath("${params.samplefile}").collect().set{ samplefile }
        poolvcfs = generatepoolvcf(chroms, samplefile)
    }

    // Create merged BAM File
    mergeBam(chroms, bamdir)
    
    // Extract SNPs from merged BAM
    snpsMergeBam(mergeBam.out.mergedbam)

    // Genotype Pool
    snpsMergeBam.out.mergevcf.join(poolvcfs)
        .map{ chrom, vcf, vcfidx, poolvcf, poolvcfidx -> tuple(chrom, vcf, poolvcf)}
        .set{ genotypes }

    genotypePool(genotypes, wcregion, bamdir)

    assignSingleCell(genotypePool.out.results.collect())
} 