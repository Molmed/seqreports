#!/bin/bash

module load bioinfo-tools
module load FastQC

PROJECT_PATH='/nobackup/private/workspace_monika/<runfolder>/Unaligned/<Project>'

SEARCH_PATTERN='*_R[1-2]_*fastq.gz'

OUTPUT_DIR='/nobackup/private/workspace_monika/Summary_reports_v2.0/fastqc_output'

find $PROJECT_PATH -name $SEARCH_PATTERN | xargs -n 1 -P 8 -I{} fastqc -o $OUTPUT_DIR {};
