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
pip install maginator
```

Furthermore, MAGinator also needs the GTDB-tk database downloaded. Here we download release 220. If you don't already have it, you can run the following:
```sh
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
tar xvzf *.tar.gz
```

## Usage

MAGinator needs 3 input files:

* The clusters.tsv files from [VAMB](https://github.com/RasmussenLab/vamb)
* A fasta file with sequences of all contigs, with unique names
* A comma-separated file giving the position of the fastq files with your sequencing reads formatted as: SampleName,PathToForwardReads,PathToReverseReads

Run MAGinator:
```sh
maginator -v vamb_clusters.tsv -r reads.csv -c contigs.fasta -o my_output -g "/path/to/GTDB-Tk/database/release220/"
```

A testset can be found in the test_data directory. 
1. Download the 3 samples used for the test at SRA: https://www.ncbi.nlm.nih.gov/sra?LinkName=bioproject_sra_all&from_uid=715601 with the ID's dfc99c_A, f9d84e_A and 221641_A
2. Change the paths to the read-files in reads.csv
3. Unzip the contigs.fasta.gz 
4. Run MAGinator

### Run on a compute cluster
MAGinator can run on compute clusters using qsub (torque), sbatch (slurm), or drmaa structures. The --cluster argument toggles the type of compute cluster infrastructure. The --cluster_info argument toggles the information given to the submission command, and it has to contain the following keywords {cores}, {memory}, {runtime}, which are used to forward resource information to the cluster.

A qsub MAGinator can for example be run with the following command (... indicates required arguments, see above):
```sh
maginator ... --cluster qsub --cluster_info "nodes=1:ppn={resources.cores},mem={resources.mem_gb}gb,walltime={resources.runtime}"
```

## Test data

A test set can be found in the maginator/test_data directory. 
1. Download the 3 samples used for the test at SRA: https://www.ncbi.nlm.nih.gov/sra?LinkName=bioproject_sra_all&from_uid=715601 with the ID's dfc99c_A, f9d84e_A and 221641_A
2. Clone repo: git clone https://github.com/Russel88/MAGinator.git
3. Change the paths to the read-files in reads.csv
4. Unzip the contigs.fasta.gz 
5. Run MAGinator

MAGinator can been run on the test data on a slurm server with the following command:
```sh
maginator --vamb_clusters clusters.tsv --reads reads.csv --contigs contigs.fasta --gtdb_db data/release220/ --output test_out --max_mem 180
```
The expected output can be found as a zipped file on Zenodo: https://doi.org/10.5281/zenodo.8279036, where MAGinator has been run on the test data (using GTDB-tk db release207_v2) on a slurm server.

On the compute cluster each job have had access to 180gb RAM, with the following time consumption: 
real	72m27.379s
user	0m18.830s
sys	1m0.454s

If you run on a smaller server you can set the parameters --max_cores and --max_mem.

## Recommended workflow 

To generate the input files to run MAGinator we have created a recommended workflow, with preprocessing, assembly and binning* of your metagenomics reads (the rules for binning have been copied from VAMB (https://github.com/RasmussenLab/vamb/blob/master/workflow/)). 
It has been setup as a snakefile in recommended_workflow/reads_to_bins.Snakefile.

The input to the workflow is the reads.csv file. The workflow can be run using snakemake:
```
snakemake --use-conda -s reads_to_bins.Snakefile --resources mem_gb=180 --config reads=reads.csv --cores 10 --printshellcmds 
```
Once the binning is done, we recommend using a tool like dRep (https://github.com/MrOlm/drep) to create the species-level clusters. The advantage of dRep is that the clustering parameters can be modified to create clusters that belong to different taxonomic levels. An R script is located in recommended workflow/MAGinator_setup.R, that will create the different input files for MAGinator adapted to the output of dRep.

Preparing data for MAGinator run
```
sed 's/@/_/g' assembly/all_assemblies.fasta > all_assemblies.fasta
sed 's/@/_/g' vamb/clusters.tsv > clusters.tsv
```

Now you are ready to run MAGinator.

## Functional Annotation

To generate the functional annotation of the genes we recommend using EggNOG mapper (https://github.com/eggnogdb/eggnog-mapper).

You can download it and try to run it on the test data
```sh
mkdir test_out/functional_annotation
emapper.py -i test/genes/all_genes_rep_seq.fasta --output test_out/functional_annotation -m diamond --cpu 38
```

The eggNOG output can be merged with clusters.tsv and further processed to obtain functional annotations of the MAG, cluster or sample levels with the following command:
```sh
(echo -e '#sample\tMAG_cluster\tMAG\tfunction'; join -1 1 -2 1 <(awk '{print $2 "\t" $1}' clusters.tsv | sort) <(tail -n +6 annotations.tsv | head -n -3 | cut -f1,15 | grep -v '\-$' | sed 's/_[[:digit:]]\+\t/\t/' | sed 's/,/\n/g' | perl -lane '{$q = $F[0] if $#F > 0; unshift(@F, $q) if $#F == 0}; print "$F[0]\t$F[1]"' | sed 's/\tko:/\t/' | sort) | awk '{print $2 "\t" $2 "\t" $3}' | sed 's/_/\t/' | sort -k1,1 -k2,2n) > MAGfunctions.tsv
```
In this case the KEGG ortholog column 15 was picked from the eggNOG-mapper output. But by cutting e.g. column number 13, one would obtain GO terms instead. Refer to the header of the eggNOG-mapper output for other available functional annotations e.g. KEGG pathways, Pfam, CAZy, COGs, etc.


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
* Map reads to the non-redundant gene catalogue 
    * Use --min_length to filter for the minimum number of basepairs that must be aligned to keep a read
    * Use --min_identity to filter for the minimum percentage of identity of mapped read to be kept
    * Use --min_map to filter for the minimum percentage of a read that has to be mapped to be kept
* Create a gene count matrix based on a signature reads approach
    * By default, MAGinator will redistribute ambiguous mapping reads based on the profile of uniquely mapping reads
    * This can be changed with the --multi option.
* Pick non-redundant genes that are only found in one MAG cluster each
* Fit signature gene model and use the resulting signature genes to get the abundance of each MAG cluster
    * Use --num_signature_genes defines the number of signature genes used for the detection of a MAG cluster
    * Use --min_mapped_signature_genes to change minimum number of signature genes to be detected in the sample to be included in the analysis
    * Use --min_samples to alter the number of samples with the MAG cluster present in order to perform signature gene refinement
* Selection of genes used for abundance calculations
    * Use --abundance_calculation to select
        * "sum": abundance is the sum of reads per bp across the total number of signature genes (num_signature_genes) used for the abundance calculation
        * "ot_trunc" - one tail truncation: abundance is the average of reads per bp across the signature genes, but excluding the most abundant signature genes as indicated by --tail-percentage 
        * "tt_trunc" - two tailed truncation: abundance is the average of reads per bp across the signature genes but excluding the most AND LEAST abundant signature genes as indicated by --tail-percentage
* Prepare for generation of phylogenies for each MAG cluster by finding outgroups and marker genes which will be used for rooting the phylogenies
* Use the read mappings to collect SNV information for each signature gene and marker gene for each sample
* Align signature and marker genes, concatenate alignments and infer phylogenetic trees for each MAG cluster
    * Use --phylo to toggle whether use fasttree (fast, approximate) or iqtree (slow, precise) to infer phylogenies
* Infer the taxonomic scope of each gene cluster. That is, at what taxonomic level are genes from a given gene cluster found in
    * Use --tax_scope_threshold to toggle the threshold for how to find the taxonomic scope consensus
* Cluster gene clusters into synteny clusters based on how often they are found adjacent on contigs


## Output

* abundance/
    * abundance_phyloseq.RData - Phyloseq object for R, with absolute abundance and taxonomic data
    * relative_abundance_phyloseq.RData - Phyloseq object for R, with relative abundance and taxonomic data
* clusters/
    * <cluster>/<bin>.fa - Fasta files with nucleotide sequence of bins
* genes/
    * all_genes.faa - Amino acid sequences of all ORFs
    * all_genes.fna - Nucletotide sequences of all ORFs
    * all_genes_nonredundant.fasta - Nucleotide sequences of gene cluster representatives
    * all_genes_cluster.tsv - Gene clusters
    * matrix/
        * gene_count_matrix.tsv - Read count for each gene cluster for each sample
        * small_gene_count_matrix.tsv - Read count matrix only containing the genes, that does not cluster across MAG cluster
    * synteny/ - Intermediate files for synteny clustering of gene clusters
* gtdbtk/
    * <cluster>/ - GTDB-tk taxonomic annotation for each VAMB cluster
* logs/ - Log files
* mapped_reads/
    * bams/ - Bam files for mapping reads to gene clusters
* phylo/
    * alignments/ - Alignments for each signature gene
    * cluster_alignments/ - Concatenated alignments for each MAG cluster
    * pileup/ - SNV information for each MAG cluster and each sample
    * trees/ - Phylogenetic trees for each MAG cluster
    * stats.tab - Mapping information such as non-N fraction, number of signature genes and marker genes, read depth, and number of bases not reaching allele frequency cutoff 
    * stats_genes.tab - Same as above but the information is split per gene
* signature_genes/ 
    * \- R data files with signature gene optimization
    * read-count_detected-genes.pdf - Figure for each MAG cluster displaying number of identified SG's in each sample along with the number of reads mapped.
* signature_reads/
   * profiles - Read count profiles with ambiguous reads redistributed based on the uniquely mapped reads profile
* tabs/
    * gene_cluster_bins.tab - Table listing which bins each gene cluster was found in
    * gene_cluster_tax_scope.tab - Table listing the taxonomic scope of each gene cluster
    * metagenomicspecies.tab - Table listing which, if any, clusters where merged in MAG cluster and the taxonomy of those
    * signature_genes_cluster.tsv - Table with the signature genes for each MAG cluster
    * synteny_clusters.tab - Table listing the synteny cluster association for the gene clusters. Gene clusters from the same synteny cluster are genomically adjacent.
    * tax_matrix.tsv - Table with taxonomy information for MAG cluster
    
