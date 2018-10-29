#! /bin/bash
#This script will run a local copy of Illuminas interop_summary to create a summary csv-file that MultiQC can read
#RunInfo.xml and runparameters.xml is needed in addition to all interop-files to run interop_summary.
#Run the script directly on the runfolder and hopefully it will run smooth.
#Note that the output file need to be extended with _interop for MultiQC to find it.

/nobackup/private/workspace_monika/interop/bin/summary \
--csv=1 180524_D00457_0254_BCC8WMANXX >> \
/nobackup/private/workspace_monika/InterOp/180524_D00457_0254_BCC8WMANXX_interop
