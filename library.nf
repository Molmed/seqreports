nextflow.preview.dsl=2

include 'processes.nf' params(params)

def get_run_folder(run_folder) {

    Channel.value(path(run_folder))
        .ifEmpty { "Error: No run folder (--run_folder) given."; exit 1 }

}

def get_project_and_reads(run_folder) {

    Channel
        .fromPath("${run_folder}/Unaligned/**.fastq.gz", maxDepth: 5 )
        .filter( ~/^.*_[^I]\d_001\.fastq\.gz$/ )
        .map {
            it.toString.startsWith('Undetermined')?
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
        multiqc_per_flowcell( path(params.run_folder).getFileName(),
            fastqc.out.map{ it[1] }.collect().ifEmpty([]),
            fastq_screen.out.map{ it[1] }.collect().ifEmpty([])),
            interop_summary.out.map{ it[1] }.collect().ifEmpty([]),
            get_QC_thresholds.out.collect(),
            get_metadata.out.collect())
        multiqc_per_project( path(params.run_folder).getFileName(),
            combine_results_by_project(fastqc.out.groupTuple(),fastq_screen.out.groupTuple()),
            get_metadata.out.collect())

}
