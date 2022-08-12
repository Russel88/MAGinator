import sys
from snakemake import shell

class Workflow(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self, snakefile, configfile):
        cmd = 'snakemake -s {snakefile} --config wd={output} --configfile {configfile} --cores {cores}'
        cmd = cmd.format(snakefile=snakefile, configfile=configfile, output=self.output, cores=self.cores)
        shell(cmd)
