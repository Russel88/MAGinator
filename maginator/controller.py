import os
import logging
import sys
import pkg_resources
import re

class Controller(object):

    def __init__(self, ap):
     
        args = ap.parse_args()

        for k,v in args.__dict__.items():
            setattr(self, k, v)
        
        # Logger
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.log_lvl)
        logging.info('Running MAGinator version {}'.format(pkg_resources.require("maginator")[0].version))

        # Force consistency
        self.output = os.path.join(self.output, '')

        # Check cluster info input
        if self.cluster is not None:
            if self.cluster_info is not None:
                if len(re.findall('{resources.cores}|{resources.mem_gb}|{resources.runtime}', self.cluster_info)) != 3 or len(re.findall('{.*?}', self.cluster_info)) != 3: 
                    logging.error('cluster_info has to contain the following special strings: {resources.cores}, {resources.mem_gb}, and {resources.runtime}')
                    sys.exit()
                else:
                    tmp_info = re.sub('{resources.cores}|{resources.mem_gb}|{resources.runtime}','',self.cluster_info)
                    if not bool(re.match('^[a-zA-Z0-9-_ =:,.]+$', tmp_info)):
                        logging.error('Invalid characters in cluster_info')
                        sys.exit()
            else:
                logging.error('cluster_info is required when running on a compute cluster')
                sys.exit()

        # Check input and output
        self.check_params()
        self.check_out()
        if not any([self.unlock, self.only_conda, self.snake is not None]):
            self.check_vamb()
            self.check_reads()
            self.check_contigs()
        self.write_params(args)

    def check_params(self):
        '''
        Return error if parameters are not allowed
        '''

        if self.tax_scope_threshold < 0.5:
            logging.error('tax_scope_threshold lower than 0.5 can lead to unexpected results')
            sys.exit()
        if self.min_af < 0.5:
            logging.error('min_af lower than 0.5 can lead to unexpected results')
            sys.exit()

    def check_out(self):
        '''
        Creates output directory unless it exists already
        '''
        try:
            os.makedirs(self.output+'logs')
        except FileExistsError:
            logging.warning('Output directory '+self.output+' already exists')

    def check_vamb(self):
        '''
        Check the vamb input file
        '''
        logging.debug('Checking VAMB input')
        with open(self.vamb_clusters, 'r') as fh:
            vamb_file = [x.strip().split() for x in fh.readlines()]
        
        # Extract info
        self.bin_set = set([x[0] for x in vamb_file])
        self.contig_set = set([x[1] for x in vamb_file])
        self.sample_set = set([re.sub('_[0-9]+$', '', x) for x in self.bin_set])
        
        try:
            self.cluster_set = set([int(re.sub('.*_', '', x)) for x in self.bin_set])
        except Exception:
            logging.error('Cluster names should be numeric. 1st column in the VAMB input file should be of the form SampleName_ClusterName')
            sys.exit()

        if not all([len(x)==2 for x in vamb_file]):
            logging.error('VAMB clusters.tsv is invalid. All lines should contain 2 tab-separated columns')
            sys.exit()            

        del self.bin_set
        del vamb_file

    def check_reads(self):
        '''
        Check that the reads file is of the correct format
        Check if read files exists
        '''
        logging.debug('Checking read info file')
        
        with open(self.reads, 'r') as fh:
            read_file = [x.strip().split(',') for x in fh.readlines()]

        if not all([len(x)==3 for x in read_file]):
            logging.error('read file is invalid. All lines should contain 3 comma-separated columns')
            sys.exit()            
        
        # Check if samples names match
        sample_set = set([x[0] for x in read_file])
        if not self.sample_set == sample_set:
            logging.error('Sample names in read file do not match those in the VAMB clusters.tsv file')
            sys.exit()
       
        # Check if read files exists
        if not all([os.path.isfile(y) for y in [x[1] for x in read_file]]):
            logging.error('Cannot find the fastq files with the forward reads')
            sys.exit()
        if not all([os.path.isfile(y) for y in [x[2] for x in read_file]]):
            logging.error('Cannot find the fastq files with the reverse reads')
            sys.exit()
            
    def check_contigs(self):
        '''
        Load contig file and check that it matches the clusters.tsv file
        '''
        logging.debug('Checking contig fasta file')
        
        contig_names = list()
        with open(self.contigs, 'r') as fh:
            for line in fh:
                if line.startswith('>'):
                    contig_names.append(line.strip().split()[0][1:])

        if not set(contig_names) == self.contig_set:
            logging.error('Contig names in fasta file do not match contig names in VAMB clusters.tsv file')
            sys.exit()

        del self.contig_set
        del contig_names

    def write_params(self, args):
        '''
        Write parameters for snakemake workflows to a file
        '''

        pars = vars(args)
        self.params = self.output+'parameters.tab'
        fh = open(self.params, 'w')
        for k, v in pars.items():
            fh.write('{}\t{}\n'.format(k, v))
        fh.close()
