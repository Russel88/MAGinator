import sys
import yaml
import re
import os
import subprocess
import logging

class Workflow(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self, snakefile, configfile):

        if self.cluster == None:
            cmd = 'snakemake --use_conda -s {snakefile} --config wd={output} --configfile {configfile} --cores {cores} --resources mem_gb={maxmemory}'
            cmd = cmd.format(snakefile=snakefile, configfile=configfile, output=self.output, cores=self.max_cores, maxmemory=self.max_mem)

        # If run on a cluster
        else:
            
            # Make logging dirs
            if self.cluster in ('qsub', 'slurm'):

                try:
                    os.mkdir(self.output + 'logs/cluster_err')
                except:
                    pass

                try:
                    os.mkdir(self.output + 'logs/cluster_out')
                except:
                    pass
           
            else:
                
                try:
                    os.mkdir(self.output + 'logs/cluster')
                except:
                    pass
                
        if self.cluster == 'qsub':
            
            cluster_cmd = 'qsub' + ' ' + self.cluster_info
            cluster_cmd = cluster_cmd + ' -e ' + self.output + 'logs/cluster_err' + ' -o ' + self.output + 'logs/cluster_out'
            cluster_cmd = '"' + cluster_cmd + '"'

            # Substitute resource information
            cluster_cmd = re.sub('{cores}', '{resources.cores}', cluster_cmd)
            cluster_cmd = re.sub('{memory}', '{resources.memory}', cluster_cmd)
            cluster_cmd = re.sub('{runtime}', '{resources.runtime}', cluster_cmd)
           
            # Format final snakemake command
            cmd = 'snakemake --use-conda --cluster {cluster_cmd} --jobs {clusterjobs} -s {snakefile} --config wd={output} --configfile {configfile} --local-cores {cores} --resources mem_gb={maxmemory}'
            cmd = cmd.format(cluster_cmd=cluster_cmd, snakefile=snakefile, configfile=configfile, output=self.output, clusterjobs=self.max_jobs, cores=self.max_cores, maxmemory=self.max_mem)

        if self.cluster == 'slurm':
            
            cluster_cmd = 'sbatch' + ' ' + self.cluster_info
            cluster_cmd = cluster_cmd + ' -e ' + self.output + 'logs/cluster_err' + ' -o ' + self.output + 'logs/cluster_out'
            cluster_cmd = "'" + cluster_cmd + "'"
            
            # Substitute resource information
            cluster_cmd = re.sub('{cores}', '{resources.cores}', cluster_cmd)
            cluster_cmd = re.sub('{memory}', '{resources.memory}', cluster_cmd)
            cluster_cmd = re.sub('{runtime}', '{resources.runtime}', cluster_cmd)

            # Format final snakemake command
            cmd = 'snakemake --use-conda --cluster {cluster_cmd} --jobs {clusterjobs} -s {snakefile} --config wd={output} --configfile {configfile} --local-cores {cores} --resources mem_gb={maxmemory}'
            cmd = cmd.format(cluster_cmd=cluster_cmd, snakefile=snakefile, configfile=configfile, output=self.output, clusterjobs=self.max_jobs, cores=self.max_cores, maxmemory=self.max_mem)
        
        if self.cluster == 'drmaa':
            
            cluster_cmd = '" ' + self.cluster_info + '"'
            
            # Substitute resource information
            cluster_cmd = re.sub('{cores}', '{resources.cores}', cluster_cmd)
            cluster_cmd = re.sub('{memory}', '{resources.memory}', cluster_cmd)
            cluster_cmd = re.sub('{runtime}', '{resources.runtime}', cluster_cmd)
            
            # Format final snakemake command
            cmd = 'snakemake --use-conda --drmaa {cluster_cmd} --jobs {clusterjobs} -s {snakefile} --config wd={output} --configfile {configfile} --drmaa-log-dir {logdir} --local-cores {cores} --resources mem_gb={maxmemory}'
            cmd = cmd.format(cluster_cmd=cluster_cmd, snakefile=snakefile, configfile=configfile, output=self.output, clusterjobs=self.max_jobs, logdir=self.output+'logs/cluster', cores=self.max_cores, maxmemory=self.max_mem)
        
        snakemake_output = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        logging.debug(snakemake_output) 

