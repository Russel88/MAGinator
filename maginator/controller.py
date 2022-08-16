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

    def check_out(self):
        '''
        Creates output directory unless it exists already
        '''
        try:
            os.mkdir(self.output)
        except FileExistsError:
            logging.warning('Output directory '+self.output+' already exists')

    def check_vamb(self):
        '''
        Check the vamb input file
        '''
        logging.debug('Checking VAMB input')
        pass


    def check_reads(self):
        '''
        Check that the reads file is of the correct format
        Check if read files exists
        '''
        logging.debug('Reading input file')
        
        with open(self.reads, 'r') as fh:
            pass
