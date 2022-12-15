import sys
import yaml
import re
import os
import subprocess
import logging
import glob

class Message(object):

    def __init__(self, obj):
        for key, val in vars(obj).items():
            setattr(self, key, val)
        
        self.commands = []

    def add(self, x):
        
        self.commands.append(x)

    def check(self):
        
        # Get latest command
        last = self.commands[-1]

        if 'workflow/filter.Snakefile' in last:
            n_clust= 0
            n_bins = 0
            for path, dirs, files in os.walk(os.path.join(self.output, 'clusters')):
                for d in dirs:
                    n_clust += 1
                for f in files:
                    n_bins += 1
            logging.info(str(n_bins) + ' bins in ' + str(n_clust) + ' VAMB clusters left after filtering')

        if 'workflow/gtdbtk.Snakefile' in last:
            n_clust = len(glob.glob(self.output+'gtdbtk/*/classify/'))
            logging.info(str(n_clust) + ' clusters could be classified')

        if 'workflow/parse_gtdbtk.Snakefile' in last:
            with open(os.path.join(self.output, 'tabs', 'metagenomicspecies.tab'), 'r') as fh:
                n_clust = 0
                n_mgs = 0
                for line in fh:
                    n_clust += len(line.strip().split('\t')[2].split(','))
                    n_mgs += 1
            logging.info(str(n_clust) + ' VAMB clusters merged into ' + str(n_mgs) + ' metagenomic species')

        if 'workflow/filter_geneclusters.Snakefile' in last:
            logging.debug('Readmapping to the clustered genes has been done.')

        if 'workflow/gene_count_mat.Snakefile' in last:
            logging.debug('A gene count matrix has been created - summarizing the readmappings for all genes in all samples')

        if 'workflow/prescreening_genes.Snakefile' in last:
            total_clust=len(glob.glob(self.output+'signature_genes/clusters/*'))
            logging.info('A total of ' + str(total_clust) + ' clusters are included in the analysis.')            

        if 'workflow/signature_genes.Snakefile' in last:
            no_SG = 0
            for f in os.listdir(os.path.join(self.output, 'signature_genes/screened/')):
                path = os.path.join(self.output, 'signature_genes/screened/', f)
                if os.path.isfile(path):
                    if (os.path.getsize(path)<60):
                        no_SG += 1
            logging.info('Too few genes was present in ' + str(no_SG) + ' MGSs in order to identify Signature Genes. The relative abundance has been set to 0.')

        if 'workflow/phylo.Snakefile' in last:
            k = 0
            for f in glob.glob(os.path.join(self.output, 'phylo', 'trees', '*.nwk')):
                if os.path.getsize(f) > 0:
                    k += 1
            logging.info('Phylogenies generated for '+str(k)+' metagenomic species')
            
