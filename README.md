# summary-report-development
In this repository we will gather useful scripts that can be used in the process of developing new sequencing report at the SNP&amp;Seq Technology platform.

This repository is not a repository for the final product and the final code. This is just a place to share scripts and ideas during the development process.

# Pre-requisites
You need to:
  - install Nextflow
  - install [Singularity (version > 2.6)](https://singularity.lbl.gov/install-linux#adding-the-mirror-and-installing)
  - Have a `.genologicsrc` file configured for the Clarity LIMS db you want to use


# Downloading fastqc-screen database
Ideally we want to do this automatically, tough right now there are problems with internet-access from the
fastq-screen containers, so this is a temporary workaround for that. Here's how you do it:

 - Download fastqc-screen from  [here](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.13.0.tar.gz).
 - Uncompress it
 - Run `fastq_screen --get_genomes`


# How to run the nextflow pipeline
Awesome, you're all set! Let's try generating reports for your favourite runfolder:
```
nextflow -c config/nextflow.config run run_multiqc_nextflow.nf \
          --runfolder ~/large_disk/180126_HSX122_0568_BHLFWLBBXX_small/ \
          --fastq_screen_db ~/large_disk/FastQ_Screen_Genomes/
```
