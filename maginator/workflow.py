import sys
import yaml
import re
import os
import subprocess
import logging
import glob

from maginator.messages import Message

class Workflow(object):

    def __init__(self, obj):
        for key, val in vars(obj).items():
            setattr(self, key, val)
        
    def add_info(self, x):
        
        # Substitute resource information
        x = re.sub('{cores}', '{resources.cores}', x)
        x = re.sub('{mem_gb}', '{resources.mem_gb}', x)
        x = re.sub('{runtime}', '{resources.runtime}', x)
   
        return x

    def run(self, snakefile):

        # Start message instance
        messages = Message(self)

        # Define core snakemake command
        cmd = ['snakemake',
               '--use-conda',
               '--latency-wait', '150',
               '-s', snakefile,
               '--config',
               'wd='+self.output,
               'reads='+self.reads,
               'contigs='+self.contigs,
               'vamb='+self.vamb_clusters,
               'params='+self.params]
        
        # If the user wants to rerun incomplete jobs
        if self.rerun_incomplete:
            cmd.append('--rerun-incomplete')
        
        # If run on server
        if self.cluster == None:
            cmd += ['--cores', str(self.max_cores)]

        # If run on a cluster
        else:
            
            # Make logging dirs
            if self.cluster in ('qsub', 'slurm'):

                try:
                    os.mkdir(self.output + 'logs/cluster_err')
                except FileExistsError:
                    pass
                try:
                    os.mkdir(self.output + 'logs/cluster_out')
                except FileExistsError:
                    pass
            else:
                try:
                    os.mkdir(self.output + 'logs/drmaa')
                except FileExistsError:
                    pass
            
            # Add cluster info to snakemake command
            cmd += ['--jobs', str(self.max_jobs),
                    '--local-cores', str(self.max_cores)]
                
        if self.cluster == 'qsub':
            # Final snakemake command
            cmd += ['--executor', 'cluster-generic', '--cluster-generic-submit-cmd', '"qsub', '-l',self.cluster_info, '-e', self.output+'logs/cluster_err', '-o', self.output+'logs/cluster_out"']
            
        if self.cluster == 'slurm':
            # Final snakemake command
            cmd += ['--executor','cluster-generic','--cluster-generic-submit-cmd', '"sbatch',self.cluster_info,'-e',self.output+'logs/cluster_err/%j.err','-o',self.output+'logs/cluster_out/%j.out"']

        if self.cluster == 'drmaa':            
            # Final snakemake command
            cmd += ['--executor cluster-generic', '--cluster-generic-submit-cmd','--drmaa-log-dir', self.output+'logs/drmaa']
        
        # Only install conda envs if only_conda
        if self.only_conda:
            logging.info('Only creating conda environments.')
            cmd.append('--conda-create-envs-only')

        # If unlocking
        if self.unlock:
            logging.info('Unlocking working directory.')
            cmd.append('--unlock')

        logging.debug(' '.join(cmd))
        # Start snakemake process and read stdout and stderr (also save in logger)
        process = subprocess.Popen(' '.join(cmd), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, text=True,shell=True)

        if self.log_lvl == 'DEBUG':
            for line in iter(process.stdout.readline, ""):
                logging.debug(line.strip())
        process.wait()
        logging.debug('Snakemake returncode: '+str(process.returncode))

        # Check returncode
        if process.returncode != 0:
            for line in iter(process.stdout.readline, ""):
                logging.error(line.strip())
            logging.error('Snakemake returncode: '+str(process.returncode))
            sys.exit()

        # Save information on what has been run
        messages.add(snakefile) 

        # If unlocking exit program
        if self.unlock:
            sys.exit()

        # Print info
        if not self.only_conda:
            messages.check()
