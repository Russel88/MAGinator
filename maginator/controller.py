import os
import logging
import sys
import pkg_resources


class Controller(object):

    def __init__(self, args):
       
        self.vamb = args.vamb
        self.reads = args.reads
        self.output = args.output
        self.keep_tmp = args.keep_tmp
        self.log_lvl = args.log_lvl
        self.cores = args.cores

        # Logger
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.log_lvl)
        logging.info('Running MAGinator version {}'.format(pkg_resources.require("maginator")[0].version))

        # Force consistency
        self.output = os.path.join(self.output, '')

        # Check input and output
        self.check_out()
        self.check_vamb()
        self.check_reads()

    def check_out(self):
        '''
        Creates output directory unless it exists already
        '''
        if False:
            try:
                os.mkdir(self.out)
            except FileExistsError:
                logging.error('Directory '+self.out+' already exists')
                sys.exit()

    def check_vamb(self):
        '''
        Check that the expected VAMB output files are present in the correct format
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
