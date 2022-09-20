[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

# MAGinator

Combining the strengths of contig and gene based methods to provide:

* Accurate abundances of species
* SNV-level resolution with phylogenetic trees on species signature genes
* Connect within-species clades with accessory genome

## Installation

```sh
conda create -n maginator -c bioconda -c conda-forge snakemake mamba
conda activate maginator
git clone https://github.com/Russel88/MAGinator.git
cd MAGinator
pip install .
```

## Usage

MAGinator needs 3 input files:

* The clusters.tsv files from [VAMB](https://github.com/RasmussenLab/vamb)
* A fasta file with sequences of all contigs, with unique names
* A comma-separated file giving the position of the read fastq file formatted as: SampleName,PathToForwardReads,PathToReverseReads

Furthermore, MAGinator also needs the GTDB-tk database version R207_v2 downloaded. If you don't already have it, you can run the following:
```sh
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz
tar xvzf gtdbtk_v2_data.tar.gz
```

Run MAGinator:
```sh
maginator --vamb_clusters clusters.tsv --reads reads.csv --contigs contigs.fasta --output my_output --gtdb_db "/path/to/GTDB-Tk/database/release207_v2/"
```

### Run on a compute cluster
MAGinator can run on compute clusters using qsub (torque), sbatch (slurm), or drmaa structures. The --cluster argument toggles the type of compute cluster infrastructure. The --cluster_info argument toggles the information given to the submission command, and it has to contain the following keywords {cores}, {memory}, {runtime}, which are used to forward resource information to the cluster.

A qsub MAGinator can for example be run with the following command (... indicates required arguments, see above):
```sh
maginator ... --cluster qsub --cluster_info "-l nodes=1:ppn={cores}:thinnode,mem={memory}gb,walltime={runtime}"
```

