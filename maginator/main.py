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
WORKFLOW_FILTER_SNAKEFILE = os.path.join(_ROOT, 'workflow', 'filter.Snakefile')
WORKFLOW_GTDBTK_SNAKEFILE = os.path.join(_ROOT, 'workflow', 'gtdbtk.Snakefile')
WORKFLOW_PARSE_GTDBTK_SNAKEFILE = os.path.join(_ROOT, 'workflow', 'parse_gtdbtk.Snakefile')

def cli():
    
    ########## Arguments ##########
    ap = argparse.ArgumentParser(description='MAGinator version {}'.format(pkg_resources.require("maginator")[0].version), add_help=False, allow_abbrev=False, formatter_class=RawTextHelpFormatter)
    
    # Required
    apr = ap.add_argument_group('required arguments')
    apr.add_argument('--vamb_clusters', help='Path to VAMB clusters.tsv file', required=True)
    apr.add_argument('--reads', help='Comma-delimited file with format: SampleName,AbsolutePathToForwardReads,AbsolutePathToReverseReads. SampleNames should match the 1st column in clusters.tsv with the pattern SampleName_{clusternumber}', required=True)
    apr.add_argument('--contigs', help='Fasta file with contig sequences. Fasta headers should match the 2nd column in the clusters.tsv file', required=True)
    apr.add_argument('--output', help='Prefix for output directory', required=True)
    apr.add_argument('--gtdb_db', help='Path to GTDB-tk database', type=str, required=True)

    # Cluster arguments
    apc = ap.add_argument_group('compute cluster arguments')
    apc.add_argument('--cluster', help='Cluster compute structure [%(default)s]', default=None, type=str, choices=[None,'qsub','slurm','drmaa'])
    apc.add_argument('--cluster_info', help='Cluster scheduler arguments when submitting cluster jobs.\nHas to contain the following special strings:\n{memory}, {cores}, and {runtime}.\nThese special strings will be substituted by maginator to indicate resources for each job.\n{memory} is substituted for the memory in GB.\n{runtime} is substituted with the time in the following format: DD:HH:MM:SS.\nCan also contain user names, groups, etc. required by the cluster scheduler', default=None, type=str)
    apc.add_argument('--max_jobs', help='Maximum number of cluster jobs [%(default)s]', default=500, type=int)
    
    # Optional
    apo = ap.add_argument_group('optional arguments')
    apo.add_argument("-h", "--help", action="help", help="show this help message and exit")
    apo.add_argument('--max_cores', help='Maximum number of cores [%(default)s]', default=40, type=int)
    apo.add_argument('--max_mem', help='Maximum memory in GB [%(default)s]', default=180, type=int)
    apo.add_argument('--log_lvl', help='Logging level [%(default)s].', default='INFO', type=str, choices=['DEBUG','INFO','WARNING','ERROR'])
    apo.add_argument('--only_conda', help='Only install conda environments, then exit', action='store_true')

    # Parameters
    app = ap.add_argument_group('parameters')
    app.add_argument('--binsize', help='Minimum bin size for inclusion [%(default)s].', default=200000, type=int)
    app.add_argument('--annotation_prevalence', help='Minimum prevalence of taxonomic assignment in a cluster of bins to call consensus [%(default)s]', default=0.75, type=float)

    ########## Workflow ##########
    master = Controller(ap)
    
    wf = Workflow(master)
    logging.info('Filtering bins')
    wf.run(snakefile=WORKFLOW_FILTER_SNAKEFILE)

    logging.info('Classifying genomes with GTDB-tk')
    wf.run(snakefile=WORKFLOW_GTDBTK_SNAKEFILE)

    logging.info('Clustering genes and parsing GTDB-tk results')
    wf.run(snakefile=WORKFLOW_PARSE_GTDBTK_SNAKEFILE)

if __name__ == '__main__':
    cli()
