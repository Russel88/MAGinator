#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pkg_resources

from os.path import join, dirname, basename
from shutil import copyfile
from argparse import RawTextHelpFormatter

from maginator.controller import Controller
from maginator.workflow import Workflow

# Snakefiles
_ROOT = os.path.abspath(os.path.dirname(__file__))
WORKFLOW_BAM_INDEX_SNAKEFILE = os.path.join(_ROOT, 'workflow', 'bam_index.Snakefile')
WORKFLOW_BAM_INDEX_CONFIG = os.path.join(_ROOT, 'workflow', 'bam_index.config.yml')

def cli():
    
    ########## Arguments ##########
    ap = argparse.ArgumentParser(description='MAGinator version {}'.format(pkg_resources.require("maginator")[0].version), add_help=False, allow_abbrev=False, formatter_class=RawTextHelpFormatter)
    
    # Required
    apr = ap.add_argument_group('required arguments')
    apr.add_argument('--vamb_clusters', help='Path to VAMB clusters.tsv file', required=True)
    apr.add_argument('--reads', help='Comma-delimited file with format: SampleName,AbsolutePathToForwardReads,AbsolutePathToReverseReads', required=True)
    apr.add_argument('--output', help='Prefix for output directory', required=True)

    # Optional
    apo = ap.add_argument_group('optional arguments')
    apo.add_argument("-h", "--help", action="help", help="show this help message and exit")
    apo.add_argument('--cluster', help='Cluster compute structure [%(default)s]', default=None, type=str, choices=[None,'qsub','slurm','drmaa'])
    apo.add_argument('--cluster_info', help='Cluster scheduler arguments when submitting cluster jobs.\nHas to contain the following special strings:\n{memory}, {cores}, and {runtime}.\nThese special strings will be substituted by maginator to indicate resources for each job.\n{memory} is substituted for the memory in GB.\n{runtime} is substituted with the time in the following format: DD:HH:MM:SS.\nCan also contain user names, groups, etc. required by the cluster scheduler', default='', type=str)
    apo.add_argument('--max_jobs', help='Maximum number of cluster jobs [%(default)s]', default=500, type=int)
    apo.add_argument('--max_cores', help='Maximum number of cores [%(default)s]', default=40, type=int)
    apo.add_argument('--max_mem', help='Maximum memory in GB [%(default)s]', default=180, type=int)
    apo.add_argument('--log_lvl', help='Logging level [%(default)s].', default='INFO', type=str, choices=['DEBUG','INFO','WARNING','ERROR'])

    ########## Workflow ##########
    master = Controller(ap)
    
    wf = Workflow(master)
    logging.info('Indexing BAM files')
    wf.run(snakefile=WORKFLOW_BAM_INDEX_SNAKEFILE, configfile = WORKFLOW_BAM_INDEX_CONFIG)

if __name__ == '__main__':
    cli()
