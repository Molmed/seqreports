//This is a simple nextflow example script that runs fastqc on a single fastq-file.
//The result is passed on to a process creating a multiqc report that finally is copied to a results directory
// The script assumes nextflow, fastqc and multiqc to be available in PATH.
//Use the following command to run the script
//nextflow run run_multiqc_nextflow.nf

//Use path to fastq-files as input
fastqs = Channel.fromPath('/summary-report-development/resources/fastqs/*.fastq.gz')

//Create directory where results should be written.
//To access and work with files, use the file method, which returns a file system object given a file path string.
results_dir = file('/summary-report-development/nextflow/results')
results_dir.mkdir()

process runFastqc {

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

    input:
    file fastqc_results

    output:
    file multiqc_results

    """
    mkdir multiqc_results
    multiqc --outdir multiqc_results fastqc_results
    """

}

//Copy results from work/ to results/
multiqc_results.subscribe { it.copyTo(results_dir) }