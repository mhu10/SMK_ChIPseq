# SMK_ChIPseq

## Workflow
![dag1](https://user-images.githubusercontent.com/38729968/233245797-ac847f5d-d4ba-4725-a80c-b7c913105b7c.svg)


## clone workflow into working directory
git clone https://github.com/mhu10/SMK_ChIPseq path/to/workdir
cd path/to/workdir

## edit config and workflow as needed
vim config.yaml or your prefered text editor


## activate snakemake
conda activate snakemake

## dry run workflow
snakemake -n

## execute workflow
snakemake --use-conda --cores 12
