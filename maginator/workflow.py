import sys
import yaml
import re

from snakemake import shell

class Workflow(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self, snakefile, configfile):

        if self.cluster == None:
            cmd = 'snakemake --use_conda -s {snakefile} --config wd={output} --configfile {configfile} --cores {cores}'
            cmd = cmd.format(snakefile=snakefile, configfile=configfile, output=self.output, cores=self.max_cores)
        
        if self.cluster == 'qsub':
            cluster_cmd = 'qsub' + ' ' + self.cluster_info
            cluster_cmd = '"' + cluster_cmd + '"'
            
            # Get resource info from config file
            with open(configfile, 'r') as fh:
                config_data = yaml.safe_load(fh)

            # Make sure cores and memory do not exceed provided maximums
            run_cores = min(config_data['cores'], self.max_cores)
            run_mem = min(config_data['memory'], self.max_mem)

            # Substitute resource information
            cluster_cmd = re.sub('{cores}', str(run_cores), cluster_cmd)
            cluster_cmd = re.sub('{memory}', str(run_mem), cluster_cmd)
            cluster_cmd = re.sub('{runtime}', config_data['runtime'], cluster_cmd)

            # Format final snakemake command
            cmd = 'snakemake --use-conda --cluster {cluster_cmd} --jobs {clusterjobs} -s {snakefile} --config wd={output} --configfile {configfile}'
            cmd = cmd.format(cluster_cmd=cluster_cmd, snakefile=snakefile, configfile=configfile, output=self.output, clusterjobs=self.max_jobs)
            
        shell(cmd)
