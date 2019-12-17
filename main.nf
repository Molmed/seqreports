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

// parameters
// params.run_folder = "/TestData/BaseSpace/180126_HSX122_0568_BHLFWLBBXX"
// runfolder = file(params.runfolder)
// runfolder_name = runfolder.getFileName()
// config_dir = file(params.config_dir)

// fastq_screen_db = file(params.fastq_screen_db)

// assets = file(params.assets)

// scripts_folder = file(params.scripts_folder)

// //Use path to fastq-file as input
// unaligned = file("$runfolder/Unaligned")
// Channel
//     .fromPath("$unaligned/**.fastq.gz", maxDepth: 3 )
//     .filter( ~/^.*_[^I]\d_001\.fastq\.gz$/)
//     .map {
//         path = file(it)
//         if (path.getFileName().toString().startsWith('Undetermined')){
//             ['NoProject', path]
//         } else {
//             // Check the folder directly under Unaligned and assume that is
//             // the project name
//             matches = path =~ /^.*\/Unaligned\/([^\/]+)\/.*\.fastq\.gz$/
//             project = matches[0][1]
//             [project, path]
//         }
//     }
//     .into{ input_fastqc; input_fastqscreen }

// include './library.nf' params(params)

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
            // path = file(it)
            // if (path.getFileName().toString().startsWith('Undetermined')){
            //     ['NoProject', path]
            // } else {
            //     // Check the folder directly under Unaligned and assume that is
            //     // the project name
            //     matches = path =~ /^.*\/Unaligned\/([^\/]+)\/.*\.fastq\.gz$/
            //     project = matches[0][1]
            //     [project, path]
            // }
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

    // publishDir file("$results_dir/flowcell_report"), mode: 'copy', overwrite: true
    // publishDir path: { additional_output_dir ? "${additional_output_dir}/flowcell_report/" : "$results_dir/flowcell_report" },
    //            saveAs: { additional_output_dir ? it : null },
    //            mode: 'copy', overwrite: true
    // errorStrategy 'ignore'

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

    // publishDir file("$results_dir/projects/"), mode: 'copy', overwrite: true
    // publishDir path: { additional_output_dir ? "${additional_output_dir}/projects/" : "$results_dir/projects" },
    //            saveAs: { additional_output_dir ? it : null },
    //            mode: 'copy', overwrite: true
    // errorStrategy 'ignore'

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

// process InteropSummary {
//
//     input:
//     file runfolder
//
//     output:
//     file runfolder_summary_interop into interop_summary_results
//
//     """
//     summary --csv=1 $runfolder > runfolder_summary_interop
//     """
// }

// process Fastqc {
//
//     input:
//     set val(project), file(fastq_file) from input_fastqc
//
//     output:
//     set val(project), file("*_fastqc.{zip,html}") into (fastqc_results_for_flowcell,fastqc_results_for_project_ungrouped)
//
//     """
//     fastqc $fastq_file
//     """
// }

// process FastqScreen {
//
//     input:
//     set val(project), file(fastq_file) from input_fastqscreen
//     file config_dir from config_dir
//     file db from fastq_screen_db
//
//     output:
//     set val(project), file("*_screen.{txt,html}") into (fastq_screen_results_for_flowcell,fastq_screen_results_for_project_ungrouped)
//
//     """
//     fastq_screen --conf $config_dir/fastq_screen.conf $fastq_file
//     """
// }

// process GetQCThresholds {
//   input:
//   file runfolder
//   file scripts_folder
//
//   output:
//   file("qc_thresholds.yaml") into qc_thresholds_result
//
//   script:
//   if (params.checkqc_config.length() > 0){
//       checkqc_config_section = "--config ${params.checkqc_config}"
//   }
//   else{
//       checkqc_config_section = ""
//   }
//
//   """
//   python $scripts_folder/get_qc_config.py --runfolder $runfolder $checkqc_config_section
//
//   """
//
// }

// process GetMetadata {
//
//     input:
//     file scripts_folder
//     file runfolder
//
//     output:
//     file 'sequencing_metadata_mqc.yaml' into sequencing_metadata_yaml
//
//     script:
//     if (params.bcl2fastq_outdir.length() > 0){
//         bcl2fastq_outdir_section = "--bcl2fastq-outdir ${params.bcl2fastq_outdir}"
//     }
//     else{
//         bcl2fastq_outdir_section = ""
//     }
//
//     """
//     python $scripts_folder/get_metadata.py --runfolder $runfolder $bcl2fastq_outdir_section &> sequencing_metadata_mqc.yaml
//     """
// }

// process MultiQCPerFlowcell {
//
//     publishDir file("$results_dir/flowcell_report"), mode: 'copy', overwrite: true
//     publishDir path: { additional_output_dir ? "${additional_output_dir}/flowcell_report/" : "$results_dir/flowcell_report" },
//                saveAs: { additional_output_dir ? it : null },
//                mode: 'copy', overwrite: true
//
//     input:
//     file (fastqc:'FastQC/*') from fastqc_results_for_flowcell.map{ it.get(1) }.collect().ifEmpty([])
//     file (fastqscreen:'FastQScreen/*') from fastq_screen_results_for_flowcell.map{ it.get(1) }.collect().ifEmpty([])
//     file (interop_summary:'Interop_summary/*') from interop_summary_results.collect().ifEmpty([])
//     file qc_thresholds from qc_thresholds_result
//     file sequencing_metadata from sequencing_metadata_yaml
//     val runfolder_name
//     file unaligned
//     file config_dir from config_dir
//     file assets from assets
//
//     output:
//     file "*multiqc_report.html" into multiqc_report_per_flowcell
//     file "*_data.zip"
//
//     """
//     multiqc \
//         --title "Flowcell report for $runfolder_name" \
//         --ignore '*/Data/Intensities/BaseCalls/L00*' \
//         --filename $runfolder_name"_multiqc_report" -z \
//         -m fastqc -m fastq_screen -m bcl2fastq -m interop -m custom_content \
//         -c $config_dir/multiqc_flowcell_config.yaml --disable_clarity -c $qc_thresholds \
//         .
//     """
//
// }

// fastqc_results_for_project_ungrouped
//     .groupTuple()
//     .map { [it.get(0), it.get(1).flatten()] }
//     .set { fastqc_results_for_project_grouped_by_project }
//
// fastq_screen_results_for_project_ungrouped
//     .groupTuple()
//     .map { [it.get(0), it.get(1).flatten()] }
//     .set { fastq_screen_results_for_project_grouped_by_project }
//
// process MultiQCPerProject {
//
//     publishDir file("$results_dir/projects/"), mode: 'copy', overwrite: true
//     publishDir path: { additional_output_dir ? "${additional_output_dir}/projects/" : "$results_dir/projects" },
//                saveAs: { additional_output_dir ? it : null },
//                mode: 'copy', overwrite: true
//
//     input:
//     set project, file(fastqc: "*") from fastqc_results_for_project_grouped_by_project
//     set project_fastq_screen, file(fastqc_screen: "*") from fastq_screen_results_for_project_grouped_by_project
//     file config_dir from config_dir
//     file sequencing_metadata from sequencing_metadata_yaml
//     val runfolder_name
//     file unaligned
//     file assets from assets
//
//     output:
//     file "$project/*multiqc_report.html" into multiqc_report_per_project
//     file "$project/*_data.zip"
//
//     """
//     multiqc \
//         --title "Report for project $project on runfolder $runfolder_name" \
//         --ignore '*/Data/Intensities/BaseCalls/L00*' \
//         --filename $project"_"$runfolder_name"_multiqc_report" -z \
//         -m fastqc -m fastq_screen -m custom_content \
//         --clarity_project $project \
//         -o $project \
//         -c $config_dir/multiqc_project_config.yaml \
//         .
//     """
//
// }
