import sys
from snakemake import shell

class Workflow(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self, snakefile, configfile):

        if self.system == 'server':
            cmd = 'snakemake -s {snakefile} --config wd={output} --configfile {configfile} --cores {cores}'
            cmd = cmd.format(snakefile=snakefile, configfile=configfile, output=self.output, cores=self.max_cores)
        
        if self.system == 'qsub':
            if len(self.cluster_info) > 0:
                cluster_cmd = 'qsub' + ' ' + self.cluster_info
            else:
                cluster_cmd = 'qsub'
            cmd = 'snakemake --cluster {cluster_cmd} --jobs {clusterjobs} -s {snakefile} --config wd={output} --configfile {configfile} --cores {cores}'
            cmd = cmd.format(cluster_cmd=cluster_cmd, snakefile=snakefile, configfile=configfile, output=self.output, cores=self.max_cores*self.max_jobs, clusterjobs=self.max_jobs)
            
        shell(cmd)
