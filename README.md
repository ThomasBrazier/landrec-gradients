# landrec-gradients


Public code for the MBE paper " Diversity in Recombination Hotspot Characteristics and Gene Structure Shape Fine-Scale Recombination Patterns in Plant Genomes" https://doi.org/10.1093/molbev/msae183



## Workflow management

This repository is managed through scripts to compute and/or copy analyses at any step and between local/remote machiens (e.g. cluster). you ust use these scripts to keep all files up to date and at the right place.

* R scripts for statistical analyses are at the root
* Reports (html/pdf) of statistical analyses are in `./Report`
* Figures are automatically generated in reports
* `./Jobs` contains bash scripts to launch actions/analyses
* `./Get` contains bash scripts to copy data between machines (local/remote)
* `./Source` contains R/Python/Bash scripts that are programs and/or functions
* `./Envs` containts all you need to install virtual env (use only `make_*.sh` scripts)


Additionnally, pipeline have been copied into this repository, such as:
* `LDhat-snakemake-pipeline`
* `Herho`


For example, to get LD data from the cluster:

```
bash Get/get_lddata.sh
```


To estimate Rho gradients for a given <dataset>:

```
Rscript Source/job_rhogradient.R <dataset>
```
or best use 
```
bash Jobs/job_rhogradient.sh
```


To assemble the complete dataset (pooled LD maps, pooled GFF, statistics):

```
bash Jobs/job_dataset_assembly.R
```

After completing the LDhat/LDhot pipeline for all species and chromosomes, you must play the jobs in the following order:

* `job_gff_rho.sh` will parse GFF files, add the estimates of LDhat/LDhot and compute various nucleotide/diversity statistics for a given species (re-reun the job on each species/dataset)
* `dataset_assembly` will compute Rho/hotspot gradients, Rho/hotspot averages and pool the complete dataset (need to run on the complete GFF rho, all species)
* `knit_finescalegradients.sh` produces the html report summarising all results
* `landscape_figures.sh` produces genomic/rho landscapes for each chromosome in the complete dataset



