#! /usr/bin/env nextflow

nextflow.preview.dsl=2
/* ####################################################

   SNP & SEQ Run folder QC pipeline

   #################################################### */

// Pipeline parameters
params.run_folder = ""
params.result_dir = "results"
params.additional_result_dir = ""
params.bcl2fastq_outdir = ""
params.checkqc_config = ""                       // See: https://github.com/Molmed/checkQC
params.assets_dir = "$baseDir/assets"
params.config_dir = "$baseDir/config"
params.script_dir = "$baseDir/bin"
params.help = false

def helpMessage() {

    log.info """

    SNP & SEQ Run folder QC pipeline.

    This workflow runs the following tools on a run folder:
        * InterOp summary    (http://illumina.github.io/interop/example_summary.html)
        * FastQC             (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * FastqScreen        (https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
        * CheckQC            (https://github.com/Molmed/checkQC)
        * MultiQC            (https://multiqc.info/)

    Usage:
        # Using parameters supplied in a config
        nextflow run -c custom.config -profile snpseq main.nf

        # Using parameters supplied on the command line
        nextflow run -profile snpseq main.nf \\
            --run_folder '/path/to/runfolder' \\
            --fastqscreen_databases '/path/to/databases' \\
            --checkqc_config '/path/to/checkqc.config' \\
            --bcl2fastq_outdir 'Unaligned'

    Mandatory parameters:
        --run_folder                        The folder to quality check

    Optional parameters:
        --result_dir                        Path to write results (default: results)
        --additional_result_dir             Additional path to write results.
        --bcl2fastq_outdir
        --checkqc_config
        --assets_dir
        --config_dir
        --script_dir

        --help                              Print this help message.

    Notes:
        * Always quote paths that are parameters to nextflow e.g. '/path/to/file'

    """
}

if (params.help || !params.run_folder){
    helpMessage()
    exit 0
}

workflow {

    main:
    get_run_folder(params.run_folder) | check_run_quality

    publish:
    check_run_quality.out.multiqc_per_project            to: "${params.result_dir}/multiqc_by_flowcell"
    check_run_quality.out.multiqc_per_flowcell           to: "${params.result_dir}/multiqc_by_project"

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the reports in your browser.\n" : "Oops .. something went wrong." )
}

def get_run_folder(run_folder) {

    Channel.value(file(run_folder))
        .ifEmpty { "Error: No run folder (--run_folder) given."; exit 1 }

}

def get_project_and_reads(run_folder) {

    Channel
        .fromPath("${run_folder}/Unaligned/**.fastq.gz", maxDepth: 5 )
        .filter( ~/^.*_[^I]\d_001\.fastq\.gz$/ )
        .map {
            it.toString.indexOf('Undetermined') > 0 ?
                ['NoProject', it] :
                [(it.toString() =~ /^.*\/Unaligned\/([^\/]+)\/.*\.fastq\.gz$/)[0][1],it]
        }

}

def combine_results_by_project (fastqc_results,fastq_screen_results) {

    fastqc_results.mix(fastq_screen_results).groupTuple().map { it -> tuple(it[0],it[1][0],it[1][1]) }

}

workflow check_run_quality {

    /* Workflow Graph

    Runfolder -> InterOp                                -> MultiQCPerFlowcell
        |     -> GetQCThresholds                        -> MultiQCPerFlowcell
        |     -> GetMetaData                            -> MultiQCPerFlowcell + MultiQCPerProject
        |
        \ -> "Undeterminded/<Project>" -> FastQC        -> MultiQCPerFlowcell + MultiQCPerProject
                                       -> FastqScreen   -> MultiQCPerFlowcell + MultiQCPerProject
    */

    get:
        run_folder

    main:
        interop_summary(run_folder)
        get_QC_thresholds(run_folder)
        get_metadata(run_folder)
        project_and_reads = get_project_and_reads(run_folder)
        fastqc(project_and_reads)
        fastq_screen(project_and_reads)
        multiqc_per_flowcell( params.run_folder,
            fastqc.out.map{ it[1] }.collect().ifEmpty([]),
            fastq_screen.out.map{ it[1] }.collect().ifEmpty([]),
            interop_summary.out.map{ it[1] }.collect().ifEmpty([]),
            get_QC_thresholds.out.collect(),
            get_metadata.out.collect())
        multiqc_per_project( params.run_folder,
            combine_results_by_project(fastqc.out.groupTuple(),fastq_screen.out.groupTuple()),
            get_metadata.out.collect())

}

// ---------------------------------------------------
// Processes
// ---------------------------------------------------

process fastqc {

    input:
    tuple project, path(fastq_file)

    output:
    tuple project, path("*_fastqc.{zip,html}")

    script:
    """
    fastqc -t ${task.cpus} $fastq_file
    """
}

process fastq_screen {

    input:
    tuple project, path(fastq_file)

    output:
    tuple project, path("*_screen.{txt,html}")

    script:
    """
    sed -E 's/^(THREADS[[:blank:]]+)[[:digit:]]+/\1${task.cpus}/' \\
        ${params.config_dir}/fastq_screen.conf > fastq_screen.conf
    fastq_screen --conf fastq_screen.conf $fastq_file
    """
}

process get_QC_thresholds {

    input:
    path runfolder

    output:
    path("qc_thresholds.yaml")

    script:
    if ( params.checkqc_config ){
        checkqc_config_section = "--config ${params.checkqc_config}"
    } else {
        checkqc_config_section = ""
    }
    """
    python ${params.script_dir}/get_qc_config.py --runfolder $runfolder \\
        $checkqc_config_section
    """

}

process get_metadata {

    input:
    path runfolder

    output:
    path 'sequencing_metadata_mqc.yaml'

    script:
    if ( params.bcl2fastq_outdir ){
        bcl2fastq_outdir_section = "--bcl2fastq-outdir ${params.bcl2fastq_outdir}"
    } else {
        bcl2fastq_outdir_section = ""
    }
    """
    python ${params.script_dir}/get_metadata.py --runfolder $runfolder \\
        $bcl2fastq_outdir_section &> sequencing_metadata_mqc.yaml
    """
}

process interop_summary {

    input:
    path runfolder

    output:
    path runfolder_summary_interop

    script:
    """
    summary --csv=1 $runfolder > runfolder_summary_interop
    """
}

process multiqc_per_flowcell {

    input:
    val runfolder_name
    path ('FastQC/*')
    path ('FastQScreen/*')
    path ('Interop_summary/*')
    path qc_thresholds
    path sequencing_metadata

    output:
    path "*multiqc_report.html"
    path "*_data.zip"

    script:
    """
    multiqc \\
        --title "Flowcell report for ${runfolder_name}" \\
        --filename ${runfolder_name}_multiqc_report -z \\
        -m fastqc -m fastq_screen -m bcl2fastq -m interop -m custom_content \\
        -c ${params.config_dir}/multiqc_flowcell_config.yaml --disable_clarity \\
        -c ${qc_thresholds} \\
        .
    """

}

process multiqc_per_project {

    input:
    val runfolder_name
    tuple project, path("FastQC/*"), path("FastqScreen/*")
    path sequencing_metadata

    output:
    path "${project}/*multiqc_report.html"
    path "${project}/*_data.zip"

    script:
    """
    multiqc \\
        --title "Report for project ${project} on runfolder ${runfolder_name}" \\
        --filename ${project}_${runfolder_name}_multiqc_report -z \\
        -m fastqc -m fastq_screen -m custom_content \\
        --clarity_project ${project} \\
        -o ${project} \\
        -c ${params.config_dir}/multiqc_project_config.yaml \\
        .
    """

}

