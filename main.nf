#! /usr/bin/env nextflow

nextflow.preview.dsl=2
/* ####################################################

   SNP & SEQ Run folder QC pipeline

   #################################################### */

// Pipeline parameters
params.run_folder = "/path/to/run_folder"
params.result_dir = "results"
fastqscreen_default_databases = "FastQ_Screen_Genomes"
params.fastqscreen_databases = fastqscreen_default_databases
params.bcl2fastq_outdir = "Unaligned"
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
            --checkqc_config '/path/to/checkqc.config'

    Mandatory parameters:
        --run_folder                        The folder to quality check

    Optional parameters:
        --result_dir                        Path to write results (default: results)
        --bcl2fastq_outdir                  Folder name to check for fastq.gz files and demultiplexing stats (default: Unaligned)
        --checkqc_config                    Configuration file for CheckQC
        --assets_dir                        Location of project assests (default: "\$baseDir/assets").
        --config_dir                        Location of project configuration files (default: "\$baseDir/config").
        --script_dir                        Location of project scripts (default: "\$baseDir/bin")

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
    Channel.fromPath(params.run_folder,checkIfExists:true)
        .ifEmpty { "Error: No run folder (--run_folder) given."; exit 1 }
        .set {run_folder}
    check_run_quality(run_folder)

    publish:
    check_run_quality.out.projectqc            to: "${params.result_dir}/multiqc_by_project", mode: 'copy'
    check_run_quality.out.flowcellqc           to: "${params.result_dir}/multiqc_by_flowcell", mode: 'copy'

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the reports in your browser.\n" : "Oops .. something went wrong." )
}

def get_project_and_reads(run_folder) {

    Channel
        .fromPath("${run_folder}/${params.bcl2fastq_outdir}/**.fastq.gz" )
        .filter( ~/.*_[^I]\d_001\.fastq\.gz$/ )
        .ifEmpty { "Error: No fastq files found under ${run_folder}/ !\n"; exit 1 }
        .map {
            it.toString().indexOf('Undetermined') > 0 ?
                ['NoProject', it] :
                [(it.toString() =~ /^.*\/${params.bcl2fastq_outdir}\/([^\/]+)\/.*\.fastq\.gz$/)[0][1],it]
        }

}

def combine_results_by_project (fastqc_results,fastq_screen_results) {

    fastqc_results.mix(fastq_screen_results).groupTuple().map { it -> tuple(it[0],it[1][0].flatten(),it[1][1].flatten()) }

}

workflow check_run_quality {

    /* Workflow Graph

    Runfolder -> InterOp                                -> MultiQCPerFlowcell
        |     -> GetQCThresholds                        -> MultiQCPerFlowcell
        |     -> GetMetaData                            -> MultiQCPerFlowcell + MultiQCPerProject
        |
        \ -> "[Undetermined|<Project>]" -> FastQC       -> MultiQCPerFlowcell + MultiQCPerProject
                                       -> FastqScreen   -> MultiQCPerFlowcell + MultiQCPerProject
    */

    take:
        run_folder

    main:
        interop_summary(run_folder)
        get_QC_thresholds(run_folder)
        get_metadata(run_folder)
        project_and_reads = get_project_and_reads(params.run_folder)
        fastqc(project_and_reads)
        fastq_screen(project_and_reads)
        multiqc_per_flowcell( params.run_folder,
            fastqc.out.map{ it[1] }.collect(),
            fastq_screen.out.map{ it[1] }.collect(),
            interop_summary.out.collect(),
            get_QC_thresholds.out.collect().ifEmpty([]),
            get_metadata.out.collect(),
            Channel.fromPath("${params.run_folder}/${params.bcl2fastq_outdir}/Stats/Stats.json").collect().ifEmpty([]),
            params.assets_dir)
        multiqc_per_project( params.run_folder,
            combine_results_by_project(fastqc.out.groupTuple(),fastq_screen.out.groupTuple()),
            get_metadata.out.collect(),
            params.assets_dir)

    emit:
        flowcellqc = multiqc_per_flowcell.out
        projectqc = multiqc_per_project.out

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
    if [ ! -e "${params.fastqscreen_databases}" ]; then
        fastq_screen --get_genomes
    elif [ "${params.fastqscreen_databases}" != "${fastqscreen_default_databases}" ]; then
        sed -i 's#${fastqscreen_default_databases}#${params.fastqscreen_databases}#' fastq_screen.conf
    fi
    fastq_screen --conf fastq_screen.conf $fastq_file
    """
}

process get_QC_thresholds {

    input:
    path runfolder

    output:
    path "qc_thresholds.yaml" optional true

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
    path 'runfolder_summary_interop'

    script:
    """
    summary --csv=1 $runfolder > runfolder_summary_interop
    """
}

process multiqc_per_flowcell {

    input:
    val runfolder_name              // Run folder name
    path ('FastQC/*')               // Fastqc logs
    path ('FastqScreen/*')          // Fastq screen logs
    path ('Interop_summary/*')      // Interop log
    path qc_thresholds              // Quality check thresholds (optional)
    path sequencing_metadata        // Sequencing meta data ( custom content data )
    path bcl2fastq_stats            // Bcl2Fastq logs
    path assets                     // Staged copy of assets folder

    output:
    tuple path("*multiqc_report.html"), path("*_data.zip")

    script:
    threshold_parameter = qc_thresholds ? "-c ${qc_thresholds}" : ""
    """
    RUNFOLDER=\$( basename ${runfolder_name} )
    multiqc \\
        --title "Flowcell report for \${RUNFOLDER}" \\
        --filename \${RUNFOLDER}_multiqc_report.html -z \\
        -c ${params.config_dir}/multiqc_flowcell_config.yaml \\
        ${threshold_parameter} \\
        .
    """

}

process multiqc_per_project {

    input:
    val runfolder_name
    tuple project, path("FastQC/*"), path("FastqScreen/*")
    path sequencing_metadata
    path assets                     // Staged copy of assets folder

    output:
    tuple path("${project}/*multiqc_report.html"), path("${project}/*_data.zip")

    script:
    """
    RUNFOLDER=\$( basename ${runfolder_name} )
    multiqc \\
        --title "Report for project ${project} on runfolder \${RUNFOLDER}" \\
        --filename \${RUNFOLDER}_${project}_multiqc_report.html -z \\
        -o ${project} \\
        -c ${params.config_dir}/multiqc_project_config.yaml \\
        .
    """

}
