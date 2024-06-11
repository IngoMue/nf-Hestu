# nf-Hestu
A nextflow workflow to estimate heterozygosity through a genotype likelihood based approach implemented through [`ANGSD`](https://github.com/ANGSD/angsd)

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/) (Tested on version >= 22.10.2)
2. Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10). By default nextflow will create a conda environment inside the work folder, this can take quite some time and ressources and is done without using a cluster's queuing system which can cause errors. To circumvent this, one option is to use [`Mamba`](https://github.com/mamba-org/mamba) instead which is much faster and efficient than conda. Alternatively, it is also possible to create the environment manually using the environment.yml file and change the `conda` flag within `profile.config` from the .yml file to the path in which the environment was created.
3. Download the pipeline, edit or create a config profile for the cluster you are using ( [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles) ) and run the workflow. If you want to use the existing `profile.config` which is written for users of Uppmax' Rackham cluster, remember to specify your naiss project ID (format: `naiss20XX-XX-XXX`) as well as the path to `nf-Hestu/environment.yml`

    ```bash
    nextflow run main.nf -profile custom --refseq '/path/to/refseq.fa' --chrlist '/path/to/chr.list' --bamlist '/path/to/bams.list' --outdir '/path/to/outdir/'
    ```
   (see below for a more detailed example)
   Other useful flags when running nextflow scripts:
   `-resume` will resume a workflow using cached intermediate files if a previous run was terminated early
   `-with-report` prints an execution report summarising input flags, directories, resource usage and information about each task
4. Once your run has completed successfully, clean up the intermediate files.

    ```bash
    nextflow clean -f -k
    ```
**Make sure** that there is **sufficient storage capacity** on the location where you will execute the pipeline from (and where the nextflow work directory will be stored) and that the scripts located under `/nf-Hestu/bin/` have permissions to be executed (e.g. `chmod 770`).

## Overview

The pipeline performs the following steps to estimate Heterozygosity:

* Determine depth-of-coverage (DoC) thresholds per individual using 1/3 mean DoC as minimum and 2 * mean DoC as maximum <sup>1</sup>
* Calculate site allele frequency likelihoods (SAF) per individual <sup>2</sup>
* Estimate the site frequency spectrum (SFS) per individual <sup>2</sup>
* Estimate and visualise heterozygosity based on the approach suggested by [`ANGSD`](https://popgen.dk/angsd/index.php/Heterozygosity) <sup>3,4</sup>

The following tools and scripts are used for the respective steps:
<sup>1</sup> R script `/bin/01_calcMinMax.R`
<sup>2</sup> [`ANGSD`](https://github.com/ANGSD/angsd) (version 0.940)
<sup>3</sup> R script `/bin/02_Hestu.R`
<sup>4</sup> R script `/bin/03_Hestu_pop.R`

## Input

The pipeline uses indexed `.bam` files and a reference sequence (`.fasta`, `.fa`, `.fna`, ...) as input.
Input files should be provided in a tab-separated file providing Individual names, the absolute path to the bam file, mean depth-of-coverage for this individual and optionally a population label. Examples for how these files should look like are found under `/nf-Hestu/examples/`. **Please note** that this pipeline assumes that the ancestral and reference sequence are the same.

Information that needs to be supplied either through flags or specified within `nextflow.config` includes:
* Absolute path to the reference sequence `refseq` (`"/path/to/refseq.fa"`)
* The path to a file containing all chromosomes/scaffolds that should be analysed (using the same naming as for the bam files), the file should contain one chromosome/scaffold per line (see `/nf-Hestu/examples/scaffs.list`) `chrlist` (`"/path/to/chr.list"`)
* The path to a tab-separated file containing a list of individuals, the path to each bam file, each individual's mean DoC and optionally a population label `bamlist` (`/path/to/bams.list`)
* Output directory `outdir` (`"/path/to/outdir/"`)
* Specify whether you want to generate plots that show heterozygosity per population (requires population information in the input file) `inclPopPlots` (boolean, default true)
* Specify whether you are using Mamba `UseMamba` (boolean, default true)

An example for using this workflow would look like this:
```bash
nextflow run main.nf -profile custom --refseq '/some/path/reference.fa' --chrlist '/some/path/chrs.list' --bamlist '/some/path/bams_pop.list' --outdir '/some/path/output/' -resume -with-report nf-Hestu_run1.html
```
