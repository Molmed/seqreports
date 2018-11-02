#!/bin/bash

module load bioinfo-tools fastq_screen

PROJECT_PATH='/nobackup/private/workspace_monika/<runfolder>/Unaligned/<Project>'

SEARCH_PATTERN='*_R[1-2]_*fastq.gz'

OUTPUT_DIR='/nobackup/private/workspace_monika/Summary_reports_v2.0/fastq_screen_output'

find $PROJECT_PATH -name $SEARCH_PATTERN | xargs -n 1 -P 8 -I{} fastq_screen --outdir $OUTPUT_DIR \
--conf /lupus/ngi/production/latest/conf/fastq_screen.irma.conf {};
