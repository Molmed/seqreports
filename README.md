# summary-report-development
In this repository we will gather useful scripts that can be used in the process of developing new sequencing report at the SNP&amp;Seq Technology platform.

This repository is not a repository for the final product and the final code. This is just a place to share scripts and ideas during the development process.


# Downloading fastqc-screen database
Ideally we want to do this automatically, tough right now there are problems with internet-access from the
fastq-screen containers, so this is a temporary workaround for that. Here's how you do it:

 - Download fastqc-screen from  [here](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.13.0.tar.gz).
 - Uncompress it
 - Run `fastq_screen --get_genomes`
