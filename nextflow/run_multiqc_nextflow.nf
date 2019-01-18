//This is a simple nextflow script that runs fastqc on a single fastq-file.
//The result is passed on to a process creating a multiqc report that finally is copied to a results directory
// The script assumes nextflow to be available in PATH and that Singularity images are defined.
//If you want to avoid entering the Singularity image as a command line parameter,
//you can define it in the Nextflow configuration file nextflow.config.
//Use the following command to run the script
//nextflow run run_nextflow_singularity.nf -c nextflow.config

//Use path to fastq-file as input
fastqs = Channel.fromPath('/summary-report-development/resources/fastqs/example_S1_L001_R1_001.fastq.gz')

//Create directory where results should be written.
results_dir = file('/summary-report-development/nextflow/results')
results_dir.mkdir()

process runFastqc {

    publishDir results_dir, mode: 'copy', overwrite: false

    input:
    file myFastq name 'example.fastq.gz' from fastqs

    output:
    file fastqc_results

    """
    mkdir fastqc_results
    fastqc -o fastqc_results $myFastq
    """
}


process runMultiQC {

    publishDir results_dir, mode: 'copy', overwrite: false

    input:
    file fastqc_results

    output:
    file multiqc_results

    """
    mkdir multiqc_results
    multiqc --outdir multiqc_results fastqc_results
    """

}
