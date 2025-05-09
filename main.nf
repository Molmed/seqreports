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
params.config_dir = "$baseDir/config/tool_config"
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
        --config_dir                        Location of tool configuration files (default: "\$baseDir/config/tool_config").
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

def combine_results_by_project (fastqc_results,fastq_screen_results,rrna_results) {
    // fastqc_results           // [Project, [fqcfiles1, fqcfiles2, fqcfiles3]]
    // fastq_screen_results     // [Project, [fqsfiles1, fqsfiles2, fqsfiles3]]
    // rrna_results             // [Project, [rrnafiles1, rrnafiles2, rrnafiles3]]

    fastqc_results.join(fastq_screen_results)
        .join(
            rrna_results.collectFile(keepHeader:true,skip:1,sort:true) { it -> ["${it[0]}_rrna_table.tsv", it[1]] }
            .map { it -> tuple((it.name - ~/_rrna_table.tsv/), [it]) })
    // [Project, [fqcfiles1,fqcfiles2,fqcfiles3],[fqsfiles1,fqsfiles2,fqsfiles3],[Project_rrna_table.tsv]]

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
        FASTQC(project_and_reads,
            params.config_dir)
        FASTQ_SCREEN(project_and_reads,
            params.config_dir,
            params.fastqscreen_databases)
        MULTIQC_PER_FLOWCELL( params.run_folder,
            FASTQC.out.map{ it[1] }.collect(),
            FASTQ_SCREEN.out.results.map{ it[1] }.collect(),
            FASTQ_SCREEN.out.tsv.map{ it[1] }.collectFile(keepHeader:true,skip:1,sort:true),
            INTEROP_SUMMARY.out.collect(),
            GET_QC_THRESHOLDS.out.collect().ifEmpty([]),
            GET_METADATA.out.collect(),
            Channel.fromPath("${params.run_folder}/${params.bcl2fastq_outdir}/Stats/Stats.json").collect().ifEmpty([]),
            params.assets_dir,
            params.config_dir)
        MULTIQC_PER_PROJECT( params.run_folder,
            combine_results_by_project(
                FASTQC.out.groupTuple(),
                FASTQ_SCREEN.out.results.groupTuple(),
                FASTQ_SCREEN.out.tsv),
            GET_METADATA.out.collect(),
            params.assets_dir,
            params.config_dir)
}

// ---------------------------------------------------
// Processes
// ---------------------------------------------------


process FASTQC {

    label 'high_memory_and_cpus'

    input:
    tuple val(project), path(fastq_file)
    path config_dir

    output:
    tuple val(project), path("*_results")

    script:
    """
    mkdir -p $fastq_file"_fastqc_results"
    fastqc -t ${task.cpus} -a "${config_dir}/adapter_list_fastqc.txt" -o $fastq_file"_fastqc_results" $fastq_file
    """
}

process FASTQ_SCREEN {

    label 'high_memory_and_cpus'

    input:
    tuple val(project), path(fastq_file)
    path config_dir
    path fastqscreen_databases

    output:
    tuple val(project), path("*_results"), emit: results
    tuple val(project), path("rrna.tsv"), emit: tsv

    script:
    outdir = fastq_file + "_fastq_screen_results"
    sample_name = (fastq_file.name =~ /^(.*_S\d+_L\d{3}_R\d+).*/)[0][1]
    """
    sed -E 's/^(THREADS[[:blank:]]+)[[:digit:]]+/\1${task.cpus}/' \\
        ${config_dir}/fastq_screen.conf > fastq_screen.conf
    if [ ! -e "${fastqscreen_databases}" ]; then
        fastq_screen --get_genomes
    elif [ "${fastqscreen_databases}" != "${fastqscreen_default_databases}" ]; then
        sed -i 's#${fastqscreen_default_databases}#${fastqscreen_databases}#' fastq_screen.conf
    fi
    mkdir -p $outdir
    fastq_screen --conf fastq_screen.conf --outdir $outdir $fastq_file
    
    # extract rRNA numbers for custom plotting with MultiQC
    printf \"Sample\\t\" > rrna.tsv
    grep -e '^Genome' -m1 -h $outdir/*_screen.txt >> rrna.tsv
    printf \"$sample_name\\t\" >> rrna.tsv
    grep -e '^rRNA' -h $outdir/*_screen.txt >> rrna.tsv
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
    interop_summary --csv=1 $runfolder > runfolder_summary_interop
    """
}

process MULTIQC_PER_FLOWCELL {

    publishDir "${params.result_dir}/flowcell_report", mode: 'copy', overwrite: true
    label 'high_memory_and_cpus'

    input:
    val runfolder_name              // Run folder name
    path ('FastQC/*')               // Fastqc logs
    path ('FastqScreen/*')          // Fastq screen logs
    path ('rRNA/rrna_table.tsv')    // Extracted rRNA values
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
    # making a separate file to use for plotting in MultiQC since custom content can only have one plot per section
    # as described here: https://multiqc.info/docs/#introduction-1
    cp rRNA/rrna_table.tsv rRNA/rrna_plot.tsv
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
    label 'high_memory_and_cpus'

    input:
    val runfolder_name
    tuple val(project), path("FastQC/*"), path("FastqScreen/*"), path("rRNA/rrna_table.tsv")
    path sequencing_metadata
    path assets                     // Staged copy of assets folder
    path config_dir                 // Staged copy of config folder

    output:
    tuple path("${project}/*multiqc_report.html"), path("${project}/*_data.zip")

    script:
    """
    # making a separate file to use for plotting in MultiQC since custom content can only have one plot per section
    # as described here: https://multiqc.info/docs/#introduction-1
    cp rRNA/rrna_table.tsv rRNA/rrna_plot.tsv
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
