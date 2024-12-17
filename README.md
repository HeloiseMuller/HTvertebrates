# HTvertebrates

Scripts used in "Phylogenetic relatedness rather than aquatic habitat fosters horizontal transfer of transposable elements in animals" by Héloïse Muller, Rosina Savisaar, Jean Peccoud, Sylvain Charlat, Clément Gilbert.

This is a branch of the scripts used in "Horizontal transfer and evolution of transposable elements in vertebrates" by Hua-Hao Zhang, Jean Peccoud, Min-Rui-Xuan Xu, Xiao-Gu Zhang, Clément Gilbert (doi: 10.1038/s41467-020-15149-4).

These scripts are publicly available to indicate how parts of the analysis were automated. There is no guaranty regarding their use. 
For those who want to use the pipeline, see below:

## Requirements
- [R](https://cran.r-project.org) 4.3+ with the following packages:
  - [data.table](https://cran.r-project.org/web/packages/data.table/) 1.16.2
  - [stringi](https://cran.r-project.org/web/packages/stringi/) 1.8.4
  - [matrixStats](https://cran.r-project.org/web/packages/matrixStats/) 1.4.1
  - [igraph](https://cran.r-project.org/web/packages/igraph/) 2.1.1
  - [ape](https://cran.r-project.org/web/packages/ape/) 5.8
  - [seqinr](https://cran.r-project.org/web/packages/seqinr/) 1.5.1
  - [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/) 1.1-3
  - [dplyr](https://cran.r-project.org/web/packages/dplyr/) 1.1.4
  - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/) 3.5.1


  (If not found, these packages are installed automatically by the pipeline.)

- [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler) 
- [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) 4.1
- [BUSCO](https://gitlab.com/ezlab/busco) 5.4
- [MMseq2] (https://github.com/soedinglab/MMseqs2) 13.4511
- [ncbi blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 2.12.0
- [diamond](https://github.com/bbuchfink/diamond) 2.0.7
- [seqtk](https://github.com/lh3/seqtk) 1.4 

The pipeline was not tested with other versions of the above programs, but more recent versions probably work.  

Hardware requirements: a linux cluster with
- ≥200 CPUs
- ≥0.5 TB of system memory
- ≥2 TB of free hard drive space 
- an internet connection to download the genome sequences (>10 MB/s is recommended). 

On this hardware, the pipeline on the 247 animal genomes should take a few months to complete.

## File description
The R scripts whose names start with numbers performed successive stages of the analysis. The purpose of each script is described by comments at the beginning of the script. 

- `HTvFunctions.R` contains functions required for the other scripts and is sourced automatically.
- `circularPlots.R` contains functions used to draw Figure 2 of the paper. 
- The remaining scripts are launched via `Rscript` (for long, CPU-intensive tasks) from the scripts whose names start with numbers. 

The files in directory `additonal_files` are required by the scripts:
- `supplementary-data1-genomes_and_accessions.txt` gives general information about the genomes and is used to download genomes sequences from ncbi.
- `ftp_links.txt` contains URL to the genome sequences, also used to download genome sequences.
- `timetree.nwk` is the timetree (newick format) used through the analysis.
- `namedClades.txt` is a table of major vertebrate clades in this tree, with their names and color codes used to make some of the paper's figures (these figures are generated by the scripts).
- `superF.txt` makes the correspondance between repeatModeler family codes (first column), TE class (2nd column) and more common TE superfamily names (3rd column). It is used in stages 15 and 16.
- `supplementary-data3-TEcomposition_per_species.txt` is generatd by the scripts and is provided with the paper, but we also provide it here if to facilitate the reproduction of the results.

The directory `demo_TeKaKs` is provided to demo the script `TEKaKs.R` (see `Demonstration of TEKaKs.R` below), but is not required to run the pipeline.


## Steps of the pipeline

SETP A : Download genomes
Script 01 is the script used by Zhang et al. (2020) to download genomes from a file containing accession links.
Script 01_2 is the script used in Muller el al. to download the selected genomes from the list of all genomes available.

STEP B : TE annotation (scripts 02 and 03) + similarity search (script 04)

STEP C : Generate dS distribution under vertical inheritance (script 5)

STEP B \& STEP C are independent; they can be run at the same time.

STEP D : Identify TE-TE hits resulting from HT
This step can be started only once both STEP A \& STEP B are done
Script 06 runs 06bis and 06tris. Similarly, script 07 runs script 07bis.
Script 07-08 is an addition that was not present in Zhang et al. (2020). This script has to be run between scripts 07 and 08. 

STEP E : Clustering
Here we used two independent methods for clustering:
	- Clustering per clade, similarly to what was done in Zhang et al. 2020.
	- Clustering per pair of species, developed for this study

Both methods have to start with script 08 to prepare the clustering. Even though one has to use the same script for both methods, one has to comment, or uncomment, the method they want to use, or not use, at the beginning of the script.
Script 08 cannot be run from bash. Starting line 177, one has to read the comments to continue.
ATTENTION Make sure you added "WAIT between each line, as explained line 181, before running the bash scripts output by script 08.

Then, there are different versions of scripts 09 and 10 depending on the method.
For clustering per clade :
	Script 09-hitClusteringRound1_perClade.R launch 09bis-iterativeFirstClustering.R
	Script 10-hitClusteringRound2_perClade.R launch 10bis
For clustering per pairs of species :
	Script 09-hitClusteringRound1_perPairs.R does not launch other scripts (step of iterativeFirstClustering are integrated in the script)
	Script 10-hitClusteringRound2_perPairs.R does not launch other scripts 

Regarding scripts 10, make sure to run clustering 10-hitClusteringRound2_perClade.R before running script 10-hitClusteringRound2_perPairs.R, as the steps in common are not repeated. Indeed,  10-hitClusteringRound2_perPairs.R does not regenerate "involvedProtOCC200dS05.self.out".

STEP F : Apply filters to keep confident hit groups only (script 11)

STEP G : analyses in common with Zhang et al. (2020)

Scripts 12 to 16 have to be run on final output of script 10-hitClusteringRound2_perClade.R, ie occ200HitGroup_perclade.txt.

STEP H : analyses specific to Muller et al.

These scripts do not follow any numbering. Those analyses have to be run on the clustering that was done independently for each pairs of species, ie on occ200HitGroup_perPairs.txt.
test_lifeStyle.R tests for an excess of transfers in the aquatic habitat.
test_phylogeneticProximity_global.R test the effect of the phylogenetic proximity globally.
test_bayesian.R runs all Bayesian analyses.


Adapting this pipeline to other datasets, hardware configuration, and automating all procedures require modifications to the code. Some parts of the analysis were not automated.

### What you need to run the pipeline

metadata.tbl is a tab delimited table that has to contain at least these three columns:
	- assembly: the name of the assembly (GCA_xxx.1); the fasta should start by this
	- species: species name written as "species genus"
	- dbBusco: the busco database to use for this genome
it can also contain any additional columns, such as taxonomical or ecological data

tree.nwk is the tree containing all the species of the dataset.
The tips are as follow: species_genus
It has to contains the divergence time, although unprecisions are no big deal.
The important is to have a correct topology and no polytomy

Scripts need those in $PATH:
seqtk, blastn, 


## Installation
In a bash-compatible terminal that can execute git, paste
```
git clone https://github.com/HeloiseMuller/HTvertebrates.git
cd HTvertebrates/
```



### Demonstration of TEKaKs.R to compute pairwise Ka and Ks on tranposable elements
We detail how to run `TEKaKs.R` on a demo dataset, but we remind that this script (as all others) is not intended for use in any other context than the study associated with the paper.

This script computes Ka, Ks, and well as overall molecular distances on pairs of homologous transposable elements (TEs), based on HSPs between these TEs and on HSPs between TEs and proteins (HSPs are not generated by this script). See the paper's method section for a description of the approach.

The hardware requirement for this demo is a Mac/Unix/Linux computer with at last 8GB of RAM, 1GB of free hard drive space, and which is able to execute R 3.4+ in a terminal. A Windows computer cannot run the script as some R functions are not supported under windows (namely those of the `parallel` R package). 

A Mac computer may have issues installing the `igraph` R package from sources, as macOS lacks a fortran compiler. However, the `igraph` package may be installed manually by specifying NOT to install packages from sources (which is not possible to do via `Rscript`).

The other programs mentioned in the `Requirements` section need not be installed for this demo.

The `demo_TeKaKs/` directory must be immediately within the `HTvertebrates/` directory. It contains the following:
- `TEhitFile.txt` is a file of TE-TE HSPs in typical blast tabular format, but only listing sequence names and HSP coordinates.
- `blastxFile.txt` is a tabular file of TE-protein HSPs. The fields indicate the TE sequence name, start and end coordinates of the HSP on this sequence (where start < end), start coordinate of the HSP on the protein and whether the TE sequence in aligned on the protein in reverse direction.
- `fastaFile.fas` is a fasta file of the TE sequences whose names are in the two aforementioned files.

The nature of these files is also explained by comments in `TEKaKs.R`.

To run the demo, paste the following in the terminal session that you used to install the pipeline:
```
Rscript TEKaKs.R demo_TeKaKs/TEhitFile.txt demo_TeKaKs/blastxFile.txt demo_TeKaKs/fastaFile.fas demo_TeKaKs/output 2
```
where `demo_TeKaKs/output` is the output folder (automatically created) and `2` is the number of CPUs to use.

The script should run in less than 5 minutes on a standard desktop PC.

Results will be found in `demo_TeKaKs/output/allKaKs.txt`. This tabular file contains the following fields:
- `hit` is an identifier for each HSP, which corresponds to the row index of each HSP in `TEhitFile.txt`.
- `ka`, `ks`, `vka` and `vks` are the results of Ka and Ks computations (see the `kaks()` function of `seqinr`), 
- `length` is the length of the alignment on which the above were computed.
- `nMut` is the number of substitutions in this alignment.
- `K80distance` and `rawDistance` are molecular distances (according to Kimura 1980 or without any correction) between sequences in the HSP. These are computed before any of the processing required for the Ka Ks computations (the removal of certain nucleotides and codons, see the method section of the paper).


## Output of the pipeline published in Zhang et al. (2020)
More than 1TB of intermediate files are generated.
The final output corresponds to results of the publication (please see the publication for their description).
- `Figure2.pdf`, `Figure3.pdf` and `Figure4.pdf` are produced at stages 14, 15 and 16 respectively. They correspond to figures of the main text
- `FigureS1.pdf` is generated at stage 5. It corresponds to supplementary figure 1.
- `FigureS2.pdf` is generated at stage 11. It corresponds to supplementary figure 2.
- `FigureS[3-6].pdf` are generated at stage 16. They corresponds to supplementary figures 3-6.
- `tableS1.txt` is generated at stage 15. It corresponds to supplementary table 1.
- `tableS2.txt` is generated at stage 16. It corresponds to supplementary table 2.
- `supplementary-data3-TEcomposition_per_species.txt` is generated at stage 2. 
- `supplementary-data4-retained_hits.txt` is generated stage 12. 

## Output of the pipeline published in Muller et al.


