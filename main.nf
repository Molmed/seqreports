#! /usr/bin/env nextflow

nextflow.enable.dsl=2
/* ####################################################

   seqreports: SNP & SEQ Run folder QC pipeline

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

    seqreports: SNP & SEQ Run folder QC pipeline.

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

def printVersion() {

    log.info "seqreports v${workflow.manifest.version}"

}

printVersion()

if (params.help || !params.run_folder){
    helpMessage()
    exit 0
}

workflow {

    main:
    Channel.fromPath(params.run_folder,checkIfExists:true)
        .ifEmpty { "Error: No run folder (--run_folder) given."; exit 1 }
        .set {run_folder}
    CHECK_RUN_QUALITY(run_folder)

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
    // fastqc_results           // [Project, [fqcfiles1, fqcfiles2, fqcfiles3]]
    // fastq_screen_results     // [Project, [fqsfiles1, fqsfiles2, fqsfiles3]]

    fastqc_results.mix(fastq_screen_results)
        .groupTuple()           // [Project, [[fqcfiles1, fqcfiles2, fqcfiles3],[fqsfiles1,fqsfiles2,fqsfiles3]] ]
        .map { it -> tuple(it[0],it[1][0].flatten(),it[1][1].flatten()) }
                                // [Project, [fqcfiles1,fqcfiles2,fqcfiles3],[fqsfiles1,fqsfiles2,fqsfiles3]]

}

workflow CHECK_RUN_QUALITY {

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
        INTEROP_SUMMARY(run_folder)
        GET_QC_THRESHOLDS(run_folder)
        GET_METADATA(run_folder)
        project_and_reads = get_project_and_reads(params.run_folder)
        FASTQC(project_and_reads)
        FASTQ_SCREEN(project_and_reads,
		     params.config_dir,
		     params.fastqscreen_databases)
        MULTIQC_PER_FLOWCELL( params.run_folder,
            FASTQC.out.map{ it[1] }.collect(),
            FASTQ_SCREEN.out.map{ it[1] }.collect(),
            INTEROP_SUMMARY.out.collect(),
            GET_QC_THRESHOLDS.out.collect().ifEmpty([]),
            GET_METADATA.out.collect(),
            Channel.fromPath("${params.run_folder}/${params.bcl2fastq_outdir}/Stats/Stats.json").collect().ifEmpty([]),
            params.assets_dir,
            params.config_dir)
        MULTIQC_PER_PROJECT( params.run_folder,
            combine_results_by_project(FASTQC.out.groupTuple(),FASTQ_SCREEN.out.groupTuple()),
            GET_METADATA.out.collect(),
            params.assets_dir,
            params.config_dir)

}

// ---------------------------------------------------
// Processes
// ---------------------------------------------------

process FASTQC {

    input:
    tuple val(project), path(fastq_file)

    output:
    tuple val(project), path("*_results")

    script:
    """
    mkdir -p $fastq_file"_fastqc_results"
    fastqc -t ${task.cpus} -o $fastq_file"_fastqc_results" $fastq_file
    """
}

process FASTQ_SCREEN {

    input:
    tuple val(project), path(fastq_file)
    path config_dir
    path fastqscreen_databases

    output:
    tuple val(project), path("*_results")

    script:
    """
    sed -E 's/^(THREADS[[:blank:]]+)[[:digit:]]+/\1${task.cpus}/' \\
        ${config_dir}/fastq_screen.conf > fastq_screen.conf
    if [ ! -e "${fastqscreen_databases}" ]; then
        fastq_screen --get_genomes
    elif [ "${fastqscreen_databases}" != "${fastqscreen_default_databases}" ]; then
        sed -i 's#${fastqscreen_default_databases}#${fastqscreen_databases}#' fastq_screen.conf
    fi
    mkdir -p $fastq_file"_fastq_screen_results"
    fastq_screen --conf fastq_screen.conf --outdir $fastq_file"_fastq_screen_results" $fastq_file
    """
}

process GET_QC_THRESHOLDS {

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

process GET_METADATA {

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

process INTEROP_SUMMARY {

    input:
    path runfolder

    output:
    path 'runfolder_summary_interop'

    script:
    """
    summary --csv=1 $runfolder > runfolder_summary_interop
    """
}

process MULTIQC_PER_FLOWCELL {

    publishDir "${params.result_dir}/flowcell_report", mode: 'copy', overwrite: true
    label 'high_memory'

    input:
    val runfolder_name              // Run folder name
    path ('FastQC/*')               // Fastqc logs
    path ('FastqScreen/*')          // Fastq screen logs
    path ('Interop_summary/*')      // Interop log
    path qc_thresholds              // Quality check thresholds (optional)
    path sequencing_metadata        // Sequencing meta data ( custom content data )
    path bcl2fastq_stats            // Bcl2Fastq logs
    path assets                     // Staged copy of assets folder
    path config_dir                 // Staged copy of config folder

    output:
    tuple path("*multiqc_report.html"), path("*_data.zip")

    script:
    threshold_parameter = qc_thresholds ? "-c ${qc_thresholds}" : ""
    """
    RUNFOLDER=\$( basename ${runfolder_name} )
    multiqc \\
        --title "Flowcell report for \${RUNFOLDER}" \\
        --filename \${RUNFOLDER}_multiqc_report.html -z \\
        -c ${config_dir}/multiqc_main_config.yaml \\
        -c ${config_dir}/multiqc_flowcell_config.yaml \\
        ${threshold_parameter} \\
        .
    """

}

process MULTIQC_PER_PROJECT {

    publishDir "${params.result_dir}/projects", mode: 'copy', overwrite: true
    label 'high_memory'

    input:
    val runfolder_name
    tuple val(project), path("FastQC/*"), path("FastqScreen/*")
    path sequencing_metadata
    path assets                     // Staged copy of assets folder
    path config_dir                 // Staged copy of config folder

    output:
    tuple path("${project}/*multiqc_report.html"), path("${project}/*_data.zip")

    script:
    """
    RUNFOLDER=\$( basename ${runfolder_name} )
    multiqc \\
        --title "Report for project ${project} on runfolder \${RUNFOLDER}" \\
        --filename \${RUNFOLDER}_${project}_multiqc_report.html -z \\
        -o ${project} \\
        -c ${config_dir}/multiqc_main_config.yaml \\
        -c ${config_dir}/multiqc_project_config.yaml \\
        .
    """

}
