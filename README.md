# SMK_ChIPseq

This workflow performs some genernal data process for ChIPseq


## Workflow
![dag1](https://user-images.githubusercontent.com/38729968/233245797-ac847f5d-d4ba-4725-a80c-b7c913105b7c.svg)


## Install Snakemake and Git

Please make sure that Snakemake and Git are correctly installed

Snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

Git: https://anaconda.org/anaconda/git



## clone workflow into working directory

```
git clone https://github.com/mhu10/SMK_ChIPseq path/to/workdir
```


## Edit config file and workfileas needed

./SMK_ChIPseq/config/'config.yaml

./SMK_ChIPseq/Snakefile

## Activate snakemake

```
conda activate snakemake
```

## Dry run workflow

```
snakemake -n
```

## Execute workflow

```
snakemake --use-conda --cores 12
```

## Build DAG workflow chart

```
snakemake --rulegraph | dot -Tsvg > dag1.svg
```
