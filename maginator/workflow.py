import sys
import yaml
import re
import os
import subprocess
import logging
import glob

class Workflow(object):

    def __init__(self, obj):
        for key, val in vars(obj).items():
            setattr(self, key, val)
        
        self.commands = []

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
            logging.info(str(n_bins) + ' bins in ' + str(n_clust) + ' clusters left after filtering')

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
            logging.info(str(n_clust) + ' clusters merged into ' + str(n_mgs) + ' metagenomic species')

    def add_info(self, x):
        
        # Substitute resource information
        x = re.sub('{cores}', '{resources.cores}', x)
        x = re.sub('{memory}', '{resources.memory}', x)
        x = re.sub('{runtime}', '{resources.runtime}', x)
   
        return x

    def run(self, snakefile):

        # Define core snakemake command
        cmd = ['snakemake',
               '--use-conda',
               '-s', snakefile,
               '--resources', 
               'mem_gb='+str(self.max_mem),
               '--config',
               'wd='+self.output,
               'reads='+self.reads,
               'contigs='+self.contigs,
               'vamb='+self.vamb_clusters,
               'params='+self.params]
        
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
            
            cluster_cmd = 'qsub' + ' ' + self.cluster_info
            cluster_cmd = cluster_cmd + ' -e ' + self.output + 'logs/cluster_err' + ' -o ' + self.output + 'logs/cluster_out'
            cluster_cmd = self.add_info(cluster_cmd)

            # Final snakemake command
            cmd += ['--cluster', cluster_cmd]

        if self.cluster == 'slurm':
            
            cluster_cmd = 'sbatch' + ' ' + self.cluster_info
            cluster_cmd = cluster_cmd + ' -e ' + self.output + 'logs/cluster_err' + ' -o ' + self.output + 'logs/cluster_out'
            cluster_cmd = self.add_info(cluster_cmd)
            
            # Final snakemake command
            cmd += ['--cluster', cluster_cmd]

        if self.cluster == 'drmaa':
            
            cluster_cmd = self.add_info(cluster_cmd)
            
            # Final snakemake command
            cmd += ['--cluster', cluster_cmd,
                    '--drmaa-log-dir', self.output+'logs/drmaa']
        
        # Only install conda envs if only_conda
        if self.only_conda:
            cmd.append(' --conda-create-envs-only')

        # Start snakemake process and read stdout and stderr (also save in logger)
        process = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, text=True)
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
        self.commands.append(snakefile) 

        # Print info
        if not self.only_conda:
            self.check()
