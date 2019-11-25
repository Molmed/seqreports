process fastqc {

    input:
    tuple project, path(fastq_file)

    output:
    tuple project, path("*_fastqc.{zip,html}")

    """
    fastqc -t ${task.cpus} $fastq_file
    """
}

process fastq_screen {

    input:
    tuple project, path(fastq_file)

    output:
    tuple project, path("*_screen.{txt,html}")

    """
    fastq_screen --conf ${params.config_dir}/fastq_screen.conf $fastq_file
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
    path scripts_folder
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

    """
    summary --csv=1 $runfolder > runfolder_summary_interop
    """
}

process multiQC_per_flowcell {

    // publishDir file("$results_dir/flowcell_report"), mode: 'copy', overwrite: true
    // publishDir path: { additional_output_dir ? "${additional_output_dir}/flowcell_report/" : "$results_dir/flowcell_report" },
    //            saveAs: { additional_output_dir ? it : null },
    //            mode: 'copy', overwrite: true
    // errorStrategy 'ignore'

    input:
    runfolder_name
    path ('FastQC/*')
    path ('FastQScreen/*')
    path ('Interop_summary/*')
    path qc_thresholds
    path sequencing_metadata

    output:
    path "*multiqc_report.html"
    path "*_data.zip"

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

process multiQC_per_project {

    // publishDir file("$results_dir/projects/"), mode: 'copy', overwrite: true
    // publishDir path: { additional_output_dir ? "${additional_output_dir}/projects/" : "$results_dir/projects" },
    //            saveAs: { additional_output_dir ? it : null },
    //            mode: 'copy', overwrite: true
    // errorStrategy 'ignore'

    input:
    runfolder_name
    tuple project, path("FastQC/*"), path("FastqScreen/*")
    path sequencing_metadata

    output:
    path "${project}/*multiqc_report.html"
    path "${project}/*_data.zip"

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
