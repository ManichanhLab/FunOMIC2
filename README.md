# FunOMIC2

## Introduction

Here we present FunOMIC2, a new version of the FunOMIC pipeline and databases, which includes major updates in both the databases and the pipeline in order to improve the analysis of mycobiome in shotgun metagenomic samples.
<!-- [brief explanantion of the changes??] -->

## Install FunOMIC2 pipeline
1. Download the FunOMIC2 pipeline from this repository
```{bash}
git clone https://github.com/ManichanhLab/FunOMIC2.git
```
2. Create a conda environment using the funomic2_env.yml file:
```{bash}
conda env create --file=funomic2_env.yml
```
Alternatively, manually install all the dependencies listed in the file.

3. Install the following R dependencies:
	- [KEGGREST](https://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html) (tested with version 1.46)
	- [RCurl](https://cran.r-project.org/package=RCurl) (tested with version 1.98)

4. Export path of the pipeline folder
```{bash}
export PATH=$PATH:/[your current directory]/FunOMIC2
```

## Running FunOMIC2

1. Activate the conda environment
```{bash}
conda activate funomic2
```

2. Download FunOMIC2 databases from the website: [website]
```{bash}
#bacterial decontamination database
wget [website]
#FunOMIC2 taxonomy database
#FunOMIC2 protein database
```
3. Decompress FunOMIC2 database files
```{bash}
tar -Jxvf [files]
```

4. Run FunOMIC2 pipeline (replace the elements in brackets [ ] with the appropiate values)

NOTE: make sure the paths provided are absolute paths, not relative paths.

```{bash}
FunOMIC2.sh -1 [reads_1.fastq.gz] -2 [reads_2.fastq.gz] \
-p [output_prefix] -o [output directory] \
-a [bacterial database directory] \
-b [FunOMIC2-T directory] \
-c [FunOMIC2-P directory] \
-t [threads]
```

It is recomended to previously run a quality control and decontamination pipeline for metagenomic sequencing data, such as KNEADDATA (<https://huttenhower.sph.harvard.edu/kneaddata/>).

## Citing FunOMIC2

If you use FunOMIC2, please cite:

<!-- add FunOMIC2 citation -->
Xie, Z. and Manichanh, C. (2022). FunOMIC: Pipeline with built-in Fungal Taxonomic and Functional Databases for Human Mycobiome Profiling. https://doi.org/10.1016/j.csbj.2022.07.010

As well as all third-party tools used by the pipeline.
