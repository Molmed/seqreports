#! /bash/bin


singularity exec multiqc-dev.simg multiqc -m interop -m fastqc -m fastq_screen -m bcl2fastq -o /nobackup/private/workspace_monika/Summary_reports_v2.0/results/multiqc_v1.7_output /nobackup/private/workspace_monika/Summary_reports_v2.0/input_data
