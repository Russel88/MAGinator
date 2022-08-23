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
        
        x = '"' + x + '"'
        
        # Substitute resource information
        x = re.sub('{cores}', '{resources.cores}', x)
        x = re.sub('{memory}', '{resources.memory}', x)
        x = re.sub('{runtime}', '{resources.runtime}', x)
   
        return x

    def run(self, snakefile):

        # If run on server
        if self.cluster == None:
            cmd = 'snakemake --use-conda -s {snakefile} --config wd={output} reads={reads} contigs={contigs} vamb={vamb} params={params} --cores {cores} --resources mem_gb={maxmemory}'
            cmd = cmd.format(snakefile=snakefile, output=self.output, cores=self.max_cores, maxmemory=self.max_mem, reads=self.reads, contigs=self.contigs, vamb=self.vamb_clusters, params=self.params)

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
                    os.mkdir(self.output + 'logs/cluster')
                except FileExistsError:
                    pass
                
        if self.cluster == 'qsub':
            
            cluster_cmd = 'qsub' + ' ' + self.cluster_info
            cluster_cmd = cluster_cmd + ' -e ' + self.output + 'logs/cluster_err' + ' -o ' + self.output + 'logs/cluster_out'
            cluster_cmd = self.add_info(cluster_cmd)

            # Format final snakemake command
            cmd = 'snakemake --use-conda --cluster {cluster_cmd} --jobs {clusterjobs} -s {snakefile} --config wd={output} reads={reads} contigs={contigs} vamb={vamb} params={params} --local-cores {cores} --resources mem_gb={maxmemory}'
            cmd = cmd.format(cluster_cmd=cluster_cmd, snakefile=snakefile, output=self.output, clusterjobs=self.max_jobs, cores=self.max_cores, maxmemory=self.max_mem, reads=self.reads, contigs=self.contigs, vamb=self.vamb_clusters, params=self.params)

        if self.cluster == 'slurm':
            
            cluster_cmd = 'sbatch' + ' ' + self.cluster_info
            cluster_cmd = cluster_cmd + ' -e ' + self.output + 'logs/cluster_err' + ' -o ' + self.output + 'logs/cluster_out'
            cluster_cmd = self.add_info(cluster_cmd)
            
            # Format final snakemake command
            cmd = 'snakemake --use-conda --cluster {cluster_cmd} --jobs {clusterjobs} -s {snakefile} --config wd={output} reads={reads} contigs={contigs} vamb={vamb} params={params} --local-cores {cores} --resources mem_gb={maxmemory}'
            cmd = cmd.format(cluster_cmd=cluster_cmd, snakefile=snakefile, output=self.output, clusterjobs=self.max_jobs, cores=self.max_cores, maxmemory=self.max_mem, reads=self.reads, contigs=self.contigs, vamb=self.vamb_clusters, params=self.params)
        
        if self.cluster == 'drmaa':
            
            cluster_cmd = self.add_info(cluster_cmd)
            
            # Format final snakemake command
            cmd = 'snakemake --use-conda --drmaa {cluster_cmd} --jobs {clusterjobs} -s {snakefile} --config wd={output} reads={reads} contigs={contigs} vamb={vamb} params={params} --drmaa-log-dir {logdir} --local-cores {cores} --resources mem_gb={maxmemory}'
            cmd = cmd.format(cluster_cmd=cluster_cmd, snakefile=snakefile, output=self.output, clusterjobs=self.max_jobs, logdir=self.output+'logs/cluster', cores=self.max_cores, maxmemory=self.max_mem, reads=self.reads, contigs=self.contigs, vamb=self.vamb_clusters, params=self.params)

        
        # Only install conda envs if only_conda
        if self.only_conda:
            cmd = cmd + ' --conda-create-envs-only'
        
        # Start snakemake process and read stdout and stderr (also save in logger)
        process = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True, text=True)
        process.wait()
        if self.log_lvl == 'DEBUG':
            for line in iter(process.stdout.readline, ""):
                logging.debug(line.strip())

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
