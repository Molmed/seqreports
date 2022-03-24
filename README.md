# seqreports: SNP&amp;Seq Run folder QC pipeline
This is a Nextflow pipeline for generating sequencing reports for the SNP&amp;Seq Technology platform, NGI Uppsala, SciLifelab Genomics.

## Pre-requisites
You need to:
  - install Nextflow (e.g. using conda `conda create -n nextflow-env nextflow` or downloading from [nextflow.io](https://www.nextflow.io/)).
  - install [Singularity (version > 2.6)](https://singularity.lbl.gov/install-linux#adding-the-mirror-and-installing).

Optional:
  - (currently mandatory: see known issues) Download the fastq-screen database by downloading fastq-screen from [here](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.13.0.tar.gz), extract the archive and then run `fastq_screen --get_genomes`.

## How to run the nextflow pipeline
Awesome, you're all set! Let's try generating reports for your favourite runfolder:
```bash
        # Using parameters supplied in a config (see below)
        nextflow run -c custom.config -profile snpseq,singularity main.nf

        # Using parameters supplied on the command line
        nextflow run -profile snpseq,singularity main.nf \
            --run_folder '/path/to/runfolder' \
            --fastqscreen_databases '/path/to/databases' \
            --checkqc_config '/path/to/checkqc.config'
```

### Available profiles

These are the primary config profiles:
- `dev`:          Run locally with low memory.
- `uppmax`:       Uppmax slurm profile for use on the cluster at Uppmax (note: The parameter `params.project` must be supplied).
- `snpseq`:       Run locally with greater memory available than `dev`.
- `singularity`:  Enables singularity and provides container URLs.
- `test`:         Run the pipeline using test data

Additional profiles:
- `debug`: prints out the `env` properties before executing processes.

### Supplying a custom config file

Custom config files can contain all command line parameters, nextflow parameters, and overriding options.

For example:
```
resume = true
params.run_folder = '/path/to/runfolder'
params.fastqscreen_databases = '/path/to/databases'
params.checkqc_config = '/path/to/checkqc.config'
workDir = '/path/to/temporary/storage/space'
```

## Development

There are two primary branches of this project:
- `master`: The stable release branch
- `dev`: The development and test branch, to which pull requests should be made.

Tests are run through GitHub Actions when pushing code to the repo. See instructions below on how to reproduce it locally.

To keep the python parts of the project nice and tidy, we enforce that code should be formatted according to [black](https://github.com/psf/black).
To re-format your code with black, simply run:
```
black .
```

### Running tests locally

Assuming you have installed all pre-requisites (except the fastq screen database: test data comes with a minimal version of it), you can run tests locally by following these steps:

```
# create virtual environment 
virtualenv -p python3.9 venv/   

# activate venv
source venv/bin/activate

# install dependencies
pip install -r requirements-dev.txt

# run tests
pytest tests/

# perform black formatter check
black --check .
```

## Known issues:

- Unable to download genome indicies using `fastq_screen --get_genomes` as wget within the container does not resolve the address correctly. Fastq Screen must be installed separately (e.g. with conda) and the genomes downloaded prior to running the workflow. The path to the databases must then be given using the `params.fastqscreen_databases` parameter.
