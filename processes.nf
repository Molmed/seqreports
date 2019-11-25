process InteropSummary {

    input:
    path runfolder

    output:
    path runfolder_summary_interop

    """
    summary --csv=1 $runfolder > runfolder_summary_interop
    """
}

process Fastqc {

    input:
    tuple project, path(fastq_file)

    output:
    tuple project, path("*_fastqc.{zip,html}")

    """
    fastqc $fastq_file
    """
}

process FastqScreen {

    input:
    tuple project, path(fastq_file)
    path config_dir
    path db

    output:
    tuple project, path("*_screen.{txt,html}")

    """
    fastq_screen --conf $config_dir/fastq_screen.conf $fastq_file
    """
}

process GetQCThresholds {

    input:
    path runfolder
    path scripts_folder

    output:
    path("qc_thresholds.yaml")

    script:
    if ( params.checkqc_config ){
        checkqc_config_section = "--config ${params.checkqc_config}"
    } else {
        checkqc_config_section = ""
    }

    """
    python $scripts_folder/get_qc_config.py --runfolder $runfolder \\
        $checkqc_config_section
    """

}

process GetMetadata {

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
    python $scripts_folder/get_metadata.py --runfolder $runfolder \\
        $bcl2fastq_outdir_section &> sequencing_metadata_mqc.yaml
    """
}

process MultiQCPerFlowcell {

    // publishDir file("$results_dir/flowcell_report"), mode: 'copy', overwrite: true
    // publishDir path: { additional_output_dir ? "${additional_output_dir}/flowcell_report/" : "$results_dir/flowcell_report" },
    //            saveAs: { additional_output_dir ? it : null },
    //            mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    input:
    path (fastqc:'FastQC/*')
    path (fastqscreen:'FastQScreen/*')
    path (interop_summary:'Interop_summary/*')
    path qc_thresholds
    path sequencing_metadata
    runfolder_name
    path unaligned
    path config_dir
    path assets

    output:
    path "*multiqc_report.html"
    path "*_data.zip"

    """
    multiqc \\
        --title "Flowcell report for $runfolder_name" \\
        --ignore '*/Data/Intensities/BaseCalls/L00*' \\
        --filename $runfolder_name"_multiqc_report" -z \\
        -m fastqc -m fastq_screen -m bcl2fastq -m interop -m custom_content \\
        -c $config_dir/multiqc_flowcell_config.yaml --disable_clarity -c $qc_thresholds \\
        .
    """

}

process MultiQCPerProject {

    // publishDir file("$results_dir/projects/"), mode: 'copy', overwrite: true
    // publishDir path: { additional_output_dir ? "${additional_output_dir}/projects/" : "$results_dir/projects" },
    //            saveAs: { additional_output_dir ? it : null },
    //            mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    input:
    tuple project, path(fastqc: "*")
    tuple project_fastq_screen, path(fastqc_screen: "*")
    path config_dir
    path sequencing_metadata
    runfolder_name
    path unaligned
    path assets

    output:
    path "$project/*multiqc_report.html"
    path "$project/*_data.zip"

    """
    multiqc \\
        --title "Report for project $project on runfolder $runfolder_name" \\
        --ignore '*/Data/Intensities/BaseCalls/L00*' \\
        --filename $project"_"$runfolder_name"_multiqc_report" -z \\
        -m fastqc -m fastq_screen -m custom_content \\
        --clarity_project $project \\
        -o $project \\
        -c $config_dir/multiqc_project_config.yaml \\
        .
    """

}
