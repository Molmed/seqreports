// This is a simple nextflow script that runs interop_summary,fastqc, and fastqscreen on a runfolder.
// The result is passed on to a process creating a multiqc report that finally is copied to a results directory
// The script assumes nextflow to be available in PATH and that Singularity images are defined in nextflow.config.
// Use the following command to run the script
// nextflow -c config/nextflow.config run run_multiqc_nextflow.nf \
//          --runfolder ~/large_disk/180126_HSX122_0568_BHLFWLBBXX_small/ \
//          --fastq_screen_db ~/large_disk/FastQ_Screen_Genomes/
//          --checkqc_config /home/monika/git_workspace/summary-report-development/checkqc_config_monika.yaml


// ----------------
// Input parameters
// ----------------

params.runfolder = "/TestData/BaseSpace/180126_HSX122_0568_BHLFWLBBXX"
runfolder = file(params.runfolder)

//Use path to fastq-file as input
unaligned = file("$runfolder/Unaligned")
Channel
    .fromPath("$unaligned/**.fastq.gz", maxDepth: 3 )
    .filter( ~/^.*_[^I]\d_001\.fastq\.gz$/)
    .map {
        path = file(it)
        if (path.getFileName().toString().startsWith('Undetermined')){
            ['NoProject', path]
        } else {
            // Check the folder directly under Unaligned and assume that is
            // the project name
            matches = path =~ /^.*\/Unaligned\/([^\/]+)\/.*\.fastq\.gz$/
            project = matches[0][1]
            [project, path]
        }
    }
    .into{ input_fastqc; input_fastqscreen }

params.fastq_screen_config = "config/fastq_screen.conf"
fastq_screen_config = file(params.fastq_screen_config)

params.multiqc_flowcell_config = "config/multiqc_flowcell_config.yaml"
multiqc_flowcell_config = file(params.multiqc_flowcell_config)

params.multiqc_project_config = "config/multiqc_project_config.yaml"
multiqc_project_config = file(params.multiqc_project_config)

fastq_screen_db = file(params.fastq_screen_db)

params.assets = "assets/"
assets = file(params.assets)

// ---------------------------------------------------
// Create directories where results should be written.
// ---------------------------------------------------
results_dir = file('results')


// ---------------------------------------------------
// Processes
// ---------------------------------------------------
process InteropSummary {

    input:
    file runfolder

    output:
    file runfolder_summary_interop into interop_summary_results

    """
    summary --csv=1 $runfolder > runfolder_summary_interop
    """
}

process Fastqc {

    input:
    set val(project), file(fastq_file) from input_fastqc

    output:
    set val(project), file("*_fastqc.{zip,html}") into fastqc_results

    """
    fastqc $fastq_file
    """
}

fastqc_results.into{ fastqc_results_for_flowcell;  fastqc_results_for_project_ungrouped }

process FastqScreen {

    input:
    set val(project), file(fastq_file) from input_fastqscreen
    file config from fastq_screen_config
    file db from fastq_screen_db

    output:
    set val(project), file("*_screen.{txt,html}") into fastq_screen_results

    """
    fastq_screen --conf $config $fastq_file
    """
}

fastq_screen_results.into{ fastq_screen_results_for_flowcell;  fastq_screen_results_for_project_ungrouped }

checkqc_config = file(params.checkqc_config)

process GetQCThresholds {
  input:
  file runfolder
  file checkqc_config

  output:
  file("qc_thresholds.yaml") into qc_thresholds_result

  """
  /summary-report-development/bin/get_qc_config.py --runfolder $runfolder --config $checkqc_config

  """

}

process MultiQCPerFlowcell {

    publishDir file("$results_dir/flowcell_report"), mode: 'copy', overwrite: true

    input:
    file (fastqc:'FastQC/*') from fastqc_results_for_flowcell.map{ it.get(1) }.collect().ifEmpty([])
    file (fastqscreen:'FastQScreen/*') from fastq_screen_results_for_flowcell.map{ it.get(1) }.collect().ifEmpty([])
    file (interop_summary:'Interop_summary/*') from interop_summary_results.collect().ifEmpty([])
    file qc_thresholds from qc_thresholds_result
    file runfolder
    file config from multiqc_flowcell_config
    file assets from assets

    output:
    file "*multiqc_report.html" into multiqc_report_per_flowcell
    file "*_data"

    """
    multiqc \
        --title "Flowcell report for ${runfolder.getFileName()}" \
        -m fastqc -m fastq_screen -m bcl2fastq -m interop -c $config \
        --disable_clarity -c $qc_thresholds \
        .
    """

}

fastqc_results_for_project_ungrouped
    .groupTuple()
    .map { [it.get(0), it.get(1).flatten()] }
    .set { fastqc_results_for_project_grouped_by_project }

fastq_screen_results_for_project_ungrouped
    .groupTuple()
    .map { [it.get(0), it.get(1).flatten()] }
    .set { fastq_screen_results_for_project_grouped_by_project }

process MultiQCPerProject {

    publishDir file("$results_dir/projects/"), mode: 'copy', overwrite: true

    input:
    set project, file(fastqc: "*") from fastqc_results_for_project_grouped_by_project
    set project_fastq_screen, file(fastqc_screen: "*") from fastq_screen_results_for_project_grouped_by_project
    file config from multiqc_project_config
    file runfolder
    file assets from assets

    output:
    file "$project/*multiqc_report.html" into multiqc_report_per_project
    file "$project/*_data"

    """
    multiqc \
        --title "Report for project $project on runfolder ${runfolder.getFileName()}" \
        -m fastqc -m fastq_screen \
        --clarity_project $project \
        -o $project \
        -c $config \
        .
    """

}
