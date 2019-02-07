//This is a simple nextflow script that runs interop_summary,fastqc, and fastqscreen on a runfolder.
//The result is passed on to a process creating a multiqc report that finally is copied to a results directory
// The script assumes nextflow to be available in PATH and that Singularity images are defined in nextflow.config.
//Use the following command to run the script
//nextflow run run_multiqc_nextflow.nf -c nextflow.config

params.runfolder = "/TestData/BaseSpace/180126_HSX122_0568_BHLFWLBBXX"
runfolder = file(params.runfolder)

//Use path to fastq-file as input
unaligned = file("$runfolder/Unaligned")
Channel.fromPath("$unaligned/**.fastq.gz", maxDepth: 3 ).into{ input_fastqc; input_fastqscreen }

params.fastq_screen_config = "config/fastq_screen.conf"
fastq_screen_config = file(params.fastq_screen_config)

fastq_screen_db = file(params.fastq_screen_db)

//Create directory where results should be written.
results_dir = file('results')
results_dir.mkdir()

multiqc_results_dir = file("$results_dir/MultiQC")
multiqc_results_dir.mkdir()

interop_summary_results_dir = file("$results_dir/Interop_summary")
interop_summary_results_dir.mkdir()


process runInteropSummary {

    publishDir interop_summary_results_dir, mode: 'symlink', overwrite: true

    input:
    file runfolder

    output:
    file runfolder_summary_interop into interop_summary_results


    """
    summary --csv=1 $runfolder > runfolder_summary_interop
    """
}

fastqc_results_dir = file("$results_dir/FastQC")
fastqc_results_dir.mkdir()

process runFastqc {

    publishDir fastqc_results_dir, mode: 'symlink', overwrite: true

    input:
    file myFastq from input_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    """
    fastqc $myFastq
    """
}



fastqscreen_results_dir = file("$results_dir/FastQScreen")
fastqscreen_results_dir.mkdir()

process runFastqScreen {

    publishDir fastqscreen_results_dir, mode: 'symlink', overwrite: true

    input:
    file myFastq from input_fastqscreen
    file config from fastq_screen_config
    file db from fastq_screen_db

    output:
    file "*_screen.{txt,html}" into fastqscreen_results

    """
    fastq_screen --conf $config $myFastq
    """
}

process runMultiQC {

    publishDir multiqc_results_dir, mode: 'symlink', overwrite: true

    input:
    file (fastqc:'FastQC/*') from fastqc_results.collect().ifEmpty([])
    file (fastqscreen:'FastQScreen/*') from fastqscreen_results.collect().ifEmpty([])
    file (interop_summary:'Interop_summary/*') from interop_summary_results.collect().ifEmpty([])
    file runfolder

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    """
    multiqc -m bcl2fastq -m interop -m fastqc -m fastq_screen .
    """

}

