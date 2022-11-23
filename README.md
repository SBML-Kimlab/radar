
# RADAR
Rapid Analysis and Detection tool of Antimicrobial-Resistance (RADAR)

Overview:
**RADAR** is a convenient and rapid pipeline for whole genome sequence (WGS) analysis, visualization and exploration. RADAR mainly consists of three infrastructures: **Annotation** process, **Local alignment** process, and **Visualization** process.
The RADAR pipeline takes a set of assembled bacterial strains as input ( e.g. NCBI RefSeq records or user's own data in **fasta** format).
The RADAR pipeline automatically performs genomic annotation on wgs data and searches for genes through a local alignment process to the selected database. After this, the visualization process of the genes detected in the wgs data is performed. For all three processes, RADAR uses three other published tools: Prodigal, Usearch, and Circos.

RADAR pipeline provides a pipeline for researchers unfamiliar with computing using cloud services. Therefore, google colabatory, google's cloud service, was selected and has similar performance to the local pipeline.

## Table of contents
  * [Overview of dependencies](#overview-of-dependencies)
  * [Pipeline overview](#pipeline-overview)
  * [Quick start and installing dependencies](#quick-start-and-installing-dependencies)
  * [Usage](#usage)
  * [Quick start with Colabatory](#quick-start-with-colabatory)


## Pipeline overview
![RADAR](/radar.png)

## Quick start and installing dependencies

#### 1. Download RADAR pipeline on Github and Installing python packages.
```
git clone https://github.com/SBL-Kimlab/radar.git
cd radar
pip install -r requirements.txt
```
#### 2. Installing dependencies 
##### Overview of dependencies:
  * Genome annotation: [Prodigal](https://github.com/hyattpd/Prodigal)
  * Local alignment tool: [Usearch](https://www.drive5.com/usearch/)
  * Ideogram visualization: [Circos](http://circos.ca/)
##### 2-1. Genome Annotation
```
apt-get install prodigal
```
##### 2-2. Visualization 
```
# Dependencies for Circos ideogram
apt-get install libgd-dev
apt-get install cpanminus
#Installation software
cpanm Clone Config::General Font::TTF::Font GD GD::Polyline Math::Bezier Math::Round Math::VecStat Params::Validate Readonly Regexp::Common SVG Set::IntSpan Statistics::Basic Text::Format
```
## Usage

In the RADAR pipeline, there are eight different modules in detail. Each process is performed according to defined modules. Users can directly use the individual modules as shown below, so all processes can be executed at once.


```
#Before executing the RADAR pipeline, it needs to declare /include/include.ipynb.

import os; import os.path as path
path_root = path.abspath( path.join( os.getcwd() ) )
path_local = path_root +  "/radar"; path_include = path_local +  "/include"
file_include = path_include +  "/include.ipynb"
%run $file_include
```

```
#RADAR pipeline excution 
os.chdir( path_local )
SPECIES = "" #Specify the SPECIES name (e.g. escherichia)
DB = "" #Specify Database name( e.g. 1. BARDS, 2.USER_DB) 
cutoff = 0.95 # cutoff setting

radar = amr( SPECIES, DB )
radar.method.prodigal( SPECIES) 
radar.method.db_statistics( SPECIES, DB )
radar.method.blast_method.udb_making( SPECIES, DB )
radar.method.blast_method.blastp_run( SPECIES, DB ) 
radar.method.blast_parse_method.blastp_parse( SPECIES, cutoff )
radar.method.blast_parse_method.blastp_merge( SPECIES, cutoff )
```
## Quick start with Colabatory

This notebook shows how to use the cloud-based RADAR pipeline to analyze WGS data quickly and easily. The overall structure of the pipeline is very similar to a RADAR pipeline that works locally.
* [(outdated)RADAR pipeline with colab](https://colab.research.google.com/drive/1VEGQ3Fo1G8atmQ2ZtYwwkh6M8qlUw2b3?usp=sharing).
* [RADAR pipeline with colab](https://colab.research.google.com/drive/1H0m9fSBU5xnPZEHam36XcAxVOwyBRf_g?usp=sharing).

## Reference
1. *Edgar, R. C. (2010). **Search and clustering orders of magnitude faster than BLAST.** _Bioinformatics_, _26_(19), 2460-2461.*
2. *Krzywinski, M., Schein, J., Birol, I., Connors, J., Gascoyne, R., Horsman, D., ... & Marra, M. A. (2009). **Circos: an information aesthetic for comparative genomics.** _Genome research_, _19_(9), 1639-1645.*
3. *Hyatt, D., Chen, G. L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). **Prodigal: prokaryotic gene recognition and translation initiation site identification.**_BMC bioinformatics_, _11_(1), 1-11.*
4. *Alcock, B. P., Raphenya, A. R., Lau, T. T., Tsang, K. K., Bouchard, M., Edalatmand, A., ... & McArthur, A. G. (2020). **CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database.** _Nucleic acids research_, _48_(D1), D517-D525.*
