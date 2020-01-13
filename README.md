# SNP&amp;Seq summary report pipeline
This is a Nextflow pipeline for generating sequencing reports for the SNP&amp;Seq Technology platform, NGI Uppsala, SciLifelab Genomics.

## Pre-requisites
You need to:
  - install Nextflow (e.g. using conda `conda create -n nextflow-env nextflow` or downloading from [nextflow.io](https://www.nextflow.io/)).
  - install [Singularity (version > 2.6)](https://singularity.lbl.gov/install-linux#adding-the-mirror-and-installing).
  - Have a `.genologicsrc` file configured for the Clarity LIMS db you want to use.

Optional:
  - Download the fastqc-screen database by downloading fastq-screen from [here](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.13.0.tar.gz), extract the archive and then run `fastq_screen --get_genomes`.

## How to run the nextflow pipeline
Awesome, you're all set! Let's try generating reports for your favourite runfolder:
```bash
        # Using parameters supplied in a config (see below)
        nextflow run -c custom.config -profile snpseq main.nf

        # Using parameters supplied on the command line
        nextflow run -profile snpseq main.nf \
            --run_folder '/path/to/runfolder' \
            --fastqscreen_databases '/path/to/databases' \
            --checkqc_config '/path/to/checkqc.config'
```

### Available profiles

There are two primary config profiles:
- `snpseq`: For running locally.
- `irma`: For running on uppmax cluster irma with slurm (note: The parameter `params.project` must be supplied).

Additional profiles:
- `dev`: supplies less memory.
- `debug`: prints out the `env` properties before executing processes.

### Supplying a custom config file

Custom config files can contain all command line parameters and anything else that needs to be overridden.

```
resume = true
params.run_folder = '/path/to/runfolder'
params.fastqscreen_databases = '/path/to/databases'
params.checkqc_config = '/path/to/checkqc.config'
workDir = /path/to/temporary/storage/space
```
