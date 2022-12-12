[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

# MAGinator

Combining the strengths of contig and gene based methods to provide:

* Accurate abundances of species using de novo signature genes
    * MAGinator uses a statistical model to find the best genes for calculating accurate abundances
* SNV-level resolution phylogenetic trees based on signature genes
    * MAGinator creates a phylogenetic tree for each species so you can associate your metadata with subspecies/strain level differences
* Connect accessory genome to the species annotation by getting a taxonomic scope for gene clusters
    * MAGinator clusters all ORFs into gene clusters and for each gene cluster you will know which taxonomic level it is specific to
* Improve your functional annotation by grouping your genes in synteny clusters based on genomic adjacency
    * MAGinator clusters gene clusters into synteny clusters - Syntenic genes are usually part of the same pathway or have similar functions 

## Installation

All you need for running MAGinator is snakemake and mamba. Other dependencies will be installed by snakemake automatically.

```sh
conda create -n maginator -c bioconda -c conda-forge snakemake mamba
conda activate maginator
git clone https://github.com/Russel88/MAGinator.git
cd MAGinator
pip install .
```

Furthermore, MAGinator also needs the GTDB-tk database version R207_v2 downloaded. If you don't already have it, you can run the following:
```sh
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz
tar xvzf gtdbtk_v2_data.tar.gz
```

## Usage

MAGinator needs 3 input files:

* The clusters.tsv files from [VAMB](https://github.com/RasmussenLab/vamb)
* A fasta file with sequences of all contigs, with unique names
* A comma-separated file giving the position of the fastq files with your sequencing reads formatted as: SampleName,PathToForwardReads,PathToReverseReads

Run MAGinator:
```sh
maginator -v vamb_clusters.tsv -r reads.csv -c contigs.fasta -o my_output -g "/path/to/GTDB-Tk/database/release207_v2/"
```

### Run on a compute cluster
MAGinator can run on compute clusters using qsub (torque), sbatch (slurm), or drmaa structures. The --cluster argument toggles the type of compute cluster infrastructure. The --cluster_info argument toggles the information given to the submission command, and it has to contain the following keywords {cores}, {memory}, {runtime}, which are used to forward resource information to the cluster.

A qsub MAGinator can for example be run with the following command (... indicates required arguments, see above):
```sh
maginator ... --cluster qsub --cluster_info "-l nodes=1:ppn={cores}:thinnode,mem={memory}gb,walltime={runtime}"
```

## MAGinator workflow

This is what MAGinator does with your input (if you want to see all parameters run maginator --help):
* Filter bins by size
    * Use --binsize to control the cutoff
* Run GTDB-tk to taxonomically annotate bins and call open reading frames (ORFs)
* Group your VAMB clusters into metagenomic species (MGS) based on the taxonomic annotation. (Unannotated VAMB clusters are kept in the pipeline, but left unchanged)
    * Use --no_mgs to disable this
    * Use --annotation_prevalence to change how prevalent an annotation has to be in a VAMB cluster to call taxonomic consensus
* Cluster your ORFs into gene clusters to get a non-redundant gene catalogue
    * Use --clustering_min_seq_id to toggle the clustering identity
    * Use --clustering_coverage to toggle the clustering coverage
    * Use --clustering_type to toggle whether to cluster on amino acid or nucleotide level
* Map reads to the non-redundant gene catalogue and create a matrix with gene counts for each sample
* Pick non-redundant genes that are only found in one MGS each
* Fit signature gene model and use the resulting signature genes to get the abundance of each MGS
* Prepare for generation of phylogenies for each MGS by finding outgroups and marker genes which will be used for rooting the phylogenies
* Use the read mappings to collect SNV information for each signature gene and marker gene for each sample
* Align signature and marker genes, concatenate alignments and infer phylogenetic trees for each MGS
    * Use --phylo to toggle whether use fasttree (fast, approximate) or iqtree (slow, precise) to infer phylogenies
* Infer the taxonomic scope of each gene cluster. That is, at what taxonomic level are genes from a given gene cluster found in
    * Use --tax_scope_threshold to toggle the threshold for how to find the taxonomic scope consensus
* Cluster gene clusters into synteny clusters based on how often they are found adjacent on contigs


## Output

* abundance/
    * abundance_phyloseq.RData - Phyloseq object for R, with abundance and taxonomic data
* clusters/
    * <cluster>/<bin>.fa - Fasta files with nucleotide sequence of bins
* genes/
    * all_genes.faa - Amino acid sequences of all ORFs
    * all_genes.fna - Nucletotide sequences of all ORFs
    * all_genes_nonredundant.fasta - Nucleotide sequences of gene cluster representatives
    * all_genes_cluster.tsv - Gene clusters
    * matrix/
        * gene_count_matrix.tsv - Read count for each gene cluster for each sample
    * synteny/ - Intermediate files for synteny clustering of gene clusters
* gtdbtk/
    * <cluster>/ - GTDB-tk taxonomic annotation for each VAMB cluster
* logs/ - Log files
* mapped_reads/
    * bams/ - Bam files for mapping reads to gene clusters
* phylo/
    * alignments/ - Alignments for each signature gene
    * cluster_alignments/ - Concatenated alignments for each MGS
    * pileup/ - SNV information for each MGS and each sample
    * trees/ - Phylogenetic trees for each MGS
    * stats.tab - Mapping information such as non-N fraction, number of signature genes and marker genes, read depth, and number of bases not reaching allele frequency cutoff 
    * stats_genes.tab - Same as above but the information is split per gene
* signature_genes/ - R data files with signature gene optimization
* tabs/
    * gene_cluster_bins.tab - Table listing which bins each gene cluster was found in
    * gene_cluster_tax_scope.tab - Table listing the taxonomic scope of each gene cluster
    * metagenomicspecies.tab - Table listing which, if any, clusters where merged in MGS and the taxonomy of those
    * signature_genes_cluster.tsv - Table with the signature genes for each MGS/cluster
    * synteny_clusters.tab - Table listing the synteny cluster association for the gene clusters. Gene clusters from the same synteny cluster are genomically adjacent.
    * tax_matrix.tsv - Table with taxonomy information for MGS
    
