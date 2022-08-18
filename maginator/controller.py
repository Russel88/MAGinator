import os
import logging
import sys
import pkg_resources
import re

class Controller(object):

    def __init__(self, ap):
     
        args = ap.parse_args()

        self.vamb_clusters = args.vamb_clusters
        self.reads = args.reads
        self.contigs = args.contigs
        self.output = args.output
        self.cluster = args.cluster
        self.max_jobs = args.max_jobs
        self.max_cores = args.max_cores
        self.max_mem = args.max_mem
        self.log_lvl = args.log_lvl
        self.cluster_info = args.cluster_info

        # Logger
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.log_lvl)
        logging.info('Running MAGinator version {}'.format(pkg_resources.require("maginator")[0].version))

        # Force consistency
        self.output = os.path.join(self.output, '')

        # Check cluster info input
        if self.cluster is not None:
            if len(self.cluster_info) > 0:
                if len(re.findall('{cores}|{memory}|{runtime}', self.cluster_info)) != 3 or len(re.findall('{.*?}', self.cluster_info)) != 3: 
                    logging.error('cluster_info has to contain the following special strings: {cores}, {memory}, and {runtime}')
                    sys.exit()
                else:
                    tmp_info = re.sub('{cores}|{memory}|{runtime}','',self.cluster_info)
                    if not bool(re.match('^[a-zA-Z0-9-_ =:,.]+$', tmp_info)):
                        logging.error('Invalid characters in cluster_info')
                        sys.exit()
            else:
                logging.error('cluster_info is required when running on a compute cluster')
                sys.exit()

        # Check input and output
        self.check_out()
        self.check_vamb()
        self.check_reads()
        self.check_contigs()

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
        except:
            logging.error('Cluster names should be numeric. 1st column in the VAMB input file should be of the form SampleName_ClusterName')
            sys.exit()

        if not all([len(x)==2 for x in vamb_file]):
            logging.error('VAMB clusters.tsv is invalid. All lines should contain 2 tab-separated columns')
            sys.exit()            

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
            logging.error('Cannot find fastq files the forward reads')
            sys.exit()
        if not all([os.path.isfile(y) for y in [x[2] for x in read_file]]):
            logging.error('Cannot find fastq files the reverse reads')
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
            

