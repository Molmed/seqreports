#! /bash/bin

module load bioinfo-tools MultiQC

multiqc -m interop -o /nobackup/private/workspace_monika/InterOp 180524_D00457_0254_BCC8WMANXX

#This command will create a MultiQC-report for the runfolder 180524_D00457_0254_BCC8WMANXX using the interop summary.
#The interop summary is in this case stored in the file 180524_D00457_0254_BCC8WMANXX_interop inside the runfolder.
#I'm not sure at this moment if this is the only file needed for creating the report.
#I guess other files are needed since I've specified the whole runfolder in the command.
