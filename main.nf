#! /usr/bin/env nextflow

nextflow.preview.dsl=2

/* ####################################################

   SNP & SEQ Run folder QC pipeline

   #################################################### */

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
        nextflow run -profile snpseq main.nf \
            --runfolder '/path/to/runfolder' \
            --fastqscreen_databases '/path/to/databases' \
            --checkqc_config '/path/to/checkqc.config' \
            --bcl2fastq_outdir 'Unaligned'

    Notes:
        * Always quote paths that are parameters to nextflow e.g. '/path/to/file'

    """
}

if (params.help){
    helpMessage()
    exit 0
}

// parameters
// params.run_folder = "/TestData/BaseSpace/180126_HSX122_0568_BHLFWLBBXX"
params.run_folder = ""
params.result_dir = "results"
params.additional_result_dir = ""
params.bcl2fastq_outdir = ""
params.checkqc_config = ""                       // See: https://github.com/Molmed/checkQC
params.assets_dir = "$baseDir/assets"
params.config_dir = "$baseDir/config"
params.script_dir = "$baseDir/bin"
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

include 'library.nf' params(params)

workflow {

    main:
    get_run_folder(params.run_folder) | check_run_quality

    publish:
    //checkRunQuality.out.interop_summary    to: "${params.results}/interop_summary"
    //checkRunQuality.out.fastqc             to: "${params.results}/fastqc"
    //checkRunQuality.out.fastqscreen        to: "${params.results}/fastqscreen"
    checkRunQuality.out.project_multiqc            to: "${params.result_dir}/multiqc_by_flowcell"
    checkRunQuality.out.flowcell_multiqc           to: "${params.result_dir}/multiqc_by_project"

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the reports in your browser.\n" : "Oops .. something went wrong." )
}

// ---------------------------------------------------
// Processes
// ---------------------------------------------------
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
