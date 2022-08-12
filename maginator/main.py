#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pkg_resources

from os.path import join, dirname, basename
from shutil import copyfile

from maginator.controller import Controller
from maginator.workflow import Workflow

_ROOT = os.path.abspath(os.path.dirname(__file__))
WORKFLOW_BAM_INDEX_SNAKEFILE = os.path.join(_ROOT, 'workflow', 'bam_index.Snakefile')
WORKFLOW_BAM_INDEX_CONFIG = os.path.join(_ROOT, 'workflow', 'bam_index.config.yml')

def cli():
    
    ########## Arguments ##########
    ap = argparse.ArgumentParser(description='MAGinator version {}'.format(pkg_resources.require("maginator")[0].version), add_help=False)
    
    # Required
    apr = ap.add_argument_group('required arguments')
    apr.add_argument('--vamb', help='VAMB output directory', required=True)
    apr.add_argument('--reads', help='Comma-delimited file with format: SampleName,AbsolutePathToForwardReads,AbsolutePathToReverseReads', required=True)
    apr.add_argument('--output', help='Prefix for output directory', required=True)

    # Optional
    apo = ap.add_argument_group('optional arguments')
    apo.add_argument("-h", "--help", action="help", help="show this help message and exit")
    apo.add_argument('--keep_tmp', help='Keep temporary files', action='store_true')
    apo.add_argument('--log_lvl', help='Logging level [%(default)s].', default='INFO', type=str, choices=['DEBUG','INFO','WARNING','ERROR'])
    apo.add_argument('--cores', help='Maximum number of cores [%(default)s]', default=4, type=int)

    ########## Workflow ##########
    master = Controller(ap.parse_args())
    
    wf = Workflow(master)
    wf.run(snakefile=WORKFLOW_BAM_INDEX_SNAKEFILE, configfile = WORKFLOW_BAM_INDEX_CONFIG)

if __name__ == '__main__':
    cli()
