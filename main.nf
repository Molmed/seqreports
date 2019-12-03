// This is a simple nextflow script that runs interop_summary,fastqc, and fastqscreen on a runfolder.
// The result is passed on to a process creating a multiqc report that finally is copied to a results directory
// The script assumes nextflow to be available in PATH and that Singularity images are defined in nextflow.config.
// Use the following command to run the script
// nextflow -c config/nextflow.config run run_multiqc_nextflow.nf \
//          --runfolder ~/large_disk/180126_HSX122_0568_BHLFWLBBXX_small/ \
//          --fastq_screen_db ~/large_disk/FastQ_Screen_Genomes/
//          --checkqc_config /home/monika/git_workspace/summary-report-development/checkqc_config_monika.yaml
//          --bcl2fastq_outdir Unaligned


// ----------------
// Input parameters
// ----------------

params.runfolder = "/TestData/BaseSpace/180126_HSX122_0568_BHLFWLBBXX"
runfolder = file(params.runfolder)
runfolder_name = runfolder.getFileName()

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

params.config_dir = "$baseDir/config/"
config_dir = file(params.config_dir)

fastq_screen_db = file(params.fastq_screen_db)

params.bcl2fastq_outdir = ""

params.checkqc_config = ""

params.assets = "$baseDir/assets/"
assets = file(params.assets)

params.scripts_folder = "$baseDir/bin/"
scripts_folder = file(params.scripts_folder)

// ---------------------------------------------------
// Create directories where results should be written.
// ---------------------------------------------------
params.output_dir = "results"
results_dir = file(params.output_dir)

params.additional_output_dir = ""
if (params.additional_output_dir.length()) {
    additional_output_dir = file(params.additional_output_dir)
}


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
    file config_dir from config_dir
    file db from fastq_screen_db

    output:
    set val(project), file("*_screen.{txt,html}") into fastq_screen_results

    """
    fastq_screen --conf $config_dir/fastq_screen.conf $fastq_file
    """
}

fastq_screen_results.into{ fastq_screen_results_for_flowcell;  fastq_screen_results_for_project_ungrouped }

process GetQCThresholds {
  input:
  file runfolder
  file scripts_folder

  output:
  file("qc_thresholds.yaml") into qc_thresholds_result

  script:
  if (params.checkqc_config.length() > 0){
      checkqc_config_section = "--config ${params.checkqc_config}"
  }
  else{
      checkqc_config_section = ""
  }

  """
  python $scripts_folder/get_qc_config.py --runfolder $runfolder $checkqc_config_section

  """

}

process GetMetadata {

    input:
    file scripts_folder
    file runfolder

    output:
    file 'sequencing_metadata_mqc.yaml' into sequencing_metadata_yaml

    script:
    if (params.bcl2fastq_outdir.length() > 0){
        bcl2fastq_outdir_section = "--bcl2fastq-outdir ${params.bcl2fastq_outdir}"
    }
    else{
        bcl2fastq_outdir_section = ""
    }

    """
    python $scripts_folder/get_metadata.py --runfolder $runfolder $bcl2fastq_outdir_section &> sequencing_metadata_mqc.yaml
    """
}

process MultiQCPerFlowcell {

    publishDir file("$results_dir/flowcell_report"), mode: 'copy', overwrite: true
    publishDir path: { additional_output_dir ? "${additional_output_dir}/flowcell_report/" : "$results_dir/flowcell_report" },
               saveAs: { additional_output_dir ? it : null },
               mode: 'copy', overwrite: true

    input:
    file (fastqc:'FastQC/*') from fastqc_results_for_flowcell.map{ it.get(1) }.collect().ifEmpty([])
    file (fastqscreen:'FastQScreen/*') from fastq_screen_results_for_flowcell.map{ it.get(1) }.collect().ifEmpty([])
    file (interop_summary:'Interop_summary/*') from interop_summary_results.collect().ifEmpty([])
    file qc_thresholds from qc_thresholds_result
    file sequencing_metadata from sequencing_metadata_yaml
    val runfolder_name
    file unaligned
    file config_dir from config_dir
    file assets from assets

    output:
    file "*multiqc_report.html" into multiqc_report_per_flowcell
    file "*_data.zip"

    """
    multiqc \
        --title "Flowcell report for $runfolder_name" \
        --ignore '*/Data/Intensities/BaseCalls/L00*' \
        --filename $runfolder_name"_multiqc_report" -z \
        -m fastqc -m fastq_screen -m bcl2fastq -m interop -m custom_content \
        -c $config_dir/multiqc_flowcell_config.yaml -c $qc_thresholds \
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

fastqc_results_for_project_grouped_by_project
    .mix(fastq_screen_results_for_project_grouped_by_project)
    .groupTuple()
    .map { it -> tuple(it[0],it[1][0],it[1][1])}
    .set { project_results }

process MultiQCPerProject {

    publishDir file("$results_dir/projects/"), mode: 'copy', overwrite: true
    publishDir path: { additional_output_dir ? "${additional_output_dir}/projects/" : "$results_dir/projects" },
               saveAs: { additional_output_dir ? it : null },
               mode: 'copy', overwrite: true

    input:
    set project, file('fastqc/*'), file('fastqscreen/*') from project_results
    file config_dir from config_dir
    file sequencing_metadata from sequencing_metadata_yaml
    val runfolder_name
    file unaligned
    file assets from assets

    output:
    file "$project/*multiqc_report.html" into multiqc_report_per_project
    file "$project/*_data.zip"

    """
    multiqc \
        --title "Report for project $project on runfolder $runfolder_name" \
        --ignore '*/Data/Intensities/BaseCalls/L00*' \
        --filename $project"_"$runfolder_name"_multiqc_report" -z \
        -m fastqc -m fastq_screen -m custom_content \
        -o $project \
        -c $config_dir/multiqc_project_config.yaml \
        .
    """

}
