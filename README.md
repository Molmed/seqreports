# summary-report-development
This is a nextflow pipeline for generating sequencing reports for the SNP&amp;Seq Technology platform, NGI Uppsala, SciLifelab Genomics.

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
nextflow -c config/nextflow.config run main.nf \
          --runfolder ~/large_disk/180126_HSX122_0568_BHLFWLBBXX_small/ \
          --fastq_screen_db ~/large_disk/FastQ_Screen_Genomes/
```

## Profiles
There are three different config profiles:
- `standard`: For running locally
- `dev`: Like standard but with less memory
- `irma`: For running on uppmax cluster irma with slurm

Usage:
```
nextflow run main.nf -profile dev <rest of the options>
```

### irma profile
When using the irma profile, use the `--project` parameter to specify which project should be accounted for the running time