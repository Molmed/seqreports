// This is a simple nextflow script that runs interop_summary,fastqc, and fastqscreen on a runfolder.
// The result is passed on to a process creating a multiqc report that finally is copied to a results directory
// The script assumes nextflow to be available in PATH and that Singularity images are defined in nextflow.config.
// Use the following command to run the script
// nextflow -c config/nextflow.config run run_multiqc_nextflow.nf \
//          --runfolder ~/large_disk/180126_HSX122_0568_BHLFWLBBXX_small/ \
//          --fastq_screen_db ~/large_disk/FastQ_Screen_Genomes/


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

fastq_screen_db = file(params.fastq_screen_db)


// ---------------------------------------------------
// Create directories where results should be written.
// ---------------------------------------------------
results_dir = file('results')
multiqc_results_dir = file("$results_dir/MultiQC")
interop_summary_results_dir = file("$results_dir/Interop_summary")
fastqscreen_results_dir = file("$results_dir/FastQScreen")
fastqc_results_dir = file("$results_dir/FastQC")


// ---------------------------------------------------
// Processes
// ---------------------------------------------------
process InteropSummary {

    // TODO Remove
    // publishDir interop_summary_results_dir, mode: 'symlink', overwrite: true

    input:
    file runfolder

    output:
    file runfolder_summary_interop into interop_summary_results

    """
    summary --csv=1 $runfolder > runfolder_summary_interop
    """
}

process Fastqc {

    // TODO Remove
    // publishDir fastqc_results_dir, mode: 'symlink', overwrite: true

    input:
    set val(project), file(fastq_file) from input_fastqc

    output:
    set val(project), file("*_fastqc.{zip,html}") into fastqc_results

    """
    fastqc $fastq_file
    """
}

process FastqScreen {

    // TODO Remove
    // publishDir fastqscreen_results_dir, mode: 'symlink', overwrite: true

    input:
    set val(project), file(fastq_file) from input_fastqscreen
    file config from fastq_screen_config
    file db from fastq_screen_db

    output:
    set val(project), file("*_screen.{txt,html}") into fastqscreen_results

    """
    fastq_screen --conf $config $fastq_file
    """
}

fastqc_results.into{ fastqc_results_for_flowcell;  fastqc_results_for_project_ungrouped }

process MultiQCPerFlowcell {

    publishDir file("$results_dir/flowcell_report"), mode: 'symlink', overwrite: true

    input:
    file (fastqc:'FastQC/*') from fastqc_results_for_flowcell.map{ it.get(1) }.collect().ifEmpty([])
//    file (fastqscreen:'FastQScreen/*') from fastqscreen_results.collect().ifEmpty([])
//    file (interop_summary:'Interop_summary/*') from interop_summary_results.collect().ifEmpty([])
    file runfolder

    output:
    file "*multiqc_report.html" into multiqc_report_per_flowcell
    file "*_data"

    """
    multiqc -m bcl2fastq -m interop -m fastqc -m fastq_screen .
    """

}

fastqc_results_for_project_ungrouped
    .groupTuple()
    .map { [it.get(0), it.get(1).flatten()] }
    .map { println "Project: " + it.get(0); it.get(1).each{ x -> println "\t$x" }; it }
    .set { fastqc_results_for_project_grouped_by_project }

process MultiQCPerProject {

    publishDir file("$results_dir/Projects/"), mode: 'symlink', overwrite: true

    input:
    set project, file("*") from fastqc_results_for_project_grouped_by_project

    output:
    file "$project/*multiqc_report.html" into multiqc_report_per_project
    file "$project/*_data"

    """
    multiqc --title $project -m bcl2fastq -m interop -m fastqc -m fastq_screen . -o $project
    """

}
