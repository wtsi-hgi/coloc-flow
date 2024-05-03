<p align="center">
  <img src="https://github.com/wtsi-hgi/coloc/blob/main/assets/images/Coloc_flow_logo.png" width="50%"/>
</p>

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/coloc** is a bioinformatics best-practice analysis pipeline for Colocalised GWAS with eQTLs. 
Here we have integrated:
1) COJO for conditioning each of the SNPs before performing colocalisation analysis using COLOC.
   
We are currently also adding:

3) eCAVIAR (https://pubmed.ncbi.nlm.nih.gov/27866706/)
4) SMR, HEIDI (https://www.nature.com/articles/ng.3538)

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/coloc/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read the eQTL and GWAS summary statistics
2. Dynamically determine Loci to colocalise with eQTL to reduce computational burden.
3. Paralelisid COJO conditioning and COLOC analysis


## Datasets should contain these columns - 
    'P' = "p_value",
    'P-value' = "p_value",

    'REF' = "reference_allele",
    'ALT' = "alternative_allele",

    'Allele1' = "effect_allele",
    'A1' = "effect_allele",
    'A_EFF' = "effect_allele",
    'EA' = "effect_allele",

    'Allele2' = "other_allele",
    'AX' = 'other_allele',
    'A_NONEFF' = 'other_allele',
    'NEA' = 'other_allele',

    'CHROM' = "chromosome",
    'CHR' = "chromosome",

    'ID' = "variant_id",
    'MarkerName' = "variant_id",
    'SNP' = "variant_id",
    'rsID' = "variant_id",

    'POS' = "base_pair_location",
    'pos' = "base_pair_location",
    'BP' = "base_pair_location",

    'BETA' = "beta",
    'Effect' = "beta",

    'SE' = "standard_error",
    'StdErr' = "standard_error",

    'Freq1' = 'eaf',
    'A1_FREQ' = 'eaf',
    'Freq_EFF' = 'eaf',
    'EAF' = 'eaf',

    'OBS_CT' = "N",
    'N' = 'N'


If eqtl summary statistics data is used these columns should be present in files:

"gene", "SNP", "distance_to_TSS", "p", "beta"
if there are no colnames then the order should be like this, otherwise the pipeline will not work. 

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/coloc -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

## Credits

nf-core/coloc was originally written by Matiss Ozols, Iaroslav Popov, Nicola Pirastu, Charles Solomon, .

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

....... Currently maintained by HGI.

