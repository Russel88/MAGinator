#!/usr/bin/env python3

import sys
import re
import os
import shutil
import copy
import pickle

import multiprocess as mp
import pandas as pd

from Bio import SeqIO 
from collections import defaultdict

'''
Workflow of this script:
1. Get lengths of genes
2. Read mpileup files
3. Convert to a dict of genes with dicts of samples inside
4. Split to dict cluster-wise
5. For each cluster:
5.1. Get signature genes and marker genes for cluster
5.2. Get marker genes for root
5.3. Collect data, such as total non-N characters, N genes mapped per sample, coverage per sample, allele frequency threshold exceeded per sample
5.4. Filter samples per cutoffs, see above
5.5. Write statistics to file
5.6. Write sequences to fasta files for alignment
'''


#######################################################################################################################
# Constants
const_min_af = float(snakemake.params.min_af)
const_min_depth = float(snakemake.params.min_depth)
threads = int(snakemake.resources.cores)
min_nonN = float(snakemake.params.min_nonN)
min_marker_genes = int(snakemake.params.min_marker_genes)
min_signature_genes = int(snakemake.params.min_signature_genes)
pileup_dir = snakemake.input.pileup
ref_fasta = snakemake.input.ref
root_fasta = snakemake.input.root
sig_gene_tab = snakemake.input.sig
mark_gene_tab = snakemake.input.mark
tax_tab = snakemake.input.tax
align_out = snakemake.output.align
stat_out = snakemake.output.stat
stat_gene_out = snakemake.output.gene 
pileup_dump = snakemake.output.dump

#######################################################################################################################
try:
    os.mkdir(align_out)
except Exception:
    pass

def main():
    '''
    The workflow
    '''

    # Get gene lengths
    gene_len_dict = read_refs(ref_fasta)
    
    # Get root marker genes
    root_markers = read_roots(root_fasta)

    # Sample names
    mp_files = os.listdir(pileup_dir)
    sample_names = [re.sub('^(.*).mp$', '\\1', x) for x in mp_files]

    # Check if pileup files have been read before
    if os.path.exists(pileup_dump):
        with open(pileup_dump, "rb") as fhi:
            pileup_dict_gene = pickle.load(fhi)
    else:
        # Pileup files are read
        pool = mp.Pool(threads)
        pileups = list(pool.imap(mpileup_read, [os.path.join(pileup_dir, x) for x in mp_files]))

        # Change from list of samples with dicts of genes to dict of genes with dict of samples
        pileup_dict_samp = {k: v for k, v in zip(sample_names, pileups)}
        pileup_dict_gene = inverse_dict(pileup_dict_samp)

        # Dump file for re-running
        with open(pileup_dump, "wb") as fho:
            pickle.dump(pileup_dict_gene, fho)

    # Gene info
    signature_genes = pd.read_csv(sig_gene_tab, sep='\t', header=None, names=('Gene','Cluster'))
    signature_genes.dropna(inplace=True)
    tax_matrix = pd.read_csv(tax_tab, sep='\t', header=None)
    tax_matrix.dropna(inplace=True)
    signature_genes = signature_genes[signature_genes['Cluster'].isin(tax_matrix[0])]
    marker_genes = pd.read_csv(mark_gene_tab, sep='\t', header=None, names=('ClusterNo','Root','Marker','ClusterGene','RootGene'))
    marker_genes['Cluster'] = ['Cluster'+str(x) for x in marker_genes['ClusterNo']]

    with open(stat_out, 'w') as fh:
        fh.write('Cluster\tSample\tNonN\tSignature_N\tSignature_Cov\tSignature_AF\tMarker_N\tMarker_Cov\tMarker_AF\n')
    with open(stat_gene_out, 'w') as fh:
        fh.write('Cluster\tSample\tGene\tNonN\tSignature_Cov\tMarker_Cov\tSignature_AF\tMarker_AF\n')
    
    for clust in set(signature_genes['Cluster']):

        print(clust)

        # Get genes for specific cluster    
        marker_gene_list = list(marker_genes[marker_genes['Cluster'] == clust]['ClusterGene'])
        signature_gene_list = list(signature_genes[signature_genes['Cluster'] == clust]['Gene'])

        # If the signature gene is also a marker gene, only keep it as a marker gene
        signature_gene_list = [x for x in signature_gene_list if x not in marker_gene_list]

        # Extact genes for specific cluster and put in new dicts
        marker_dict = {x: pileup_dict_gene.get(x) for x in marker_gene_list if pileup_dict_gene.get(x) is not None}
        signature_dict = {x: pileup_dict_gene.get(x) for x in signature_gene_list if pileup_dict_gene.get(x) is not None}
   
        # Get marker genes for root
        root_gene_list = list(marker_genes[marker_genes['Cluster'] == clust]['RootGene'])
        root_gene_seqs = [root_markers[x] for x in root_gene_list]

        # Remove the few genes not found in pileup (bam)
        marker_gene_list = [x for x in marker_gene_list if x in marker_dict.keys()]
        signature_gene_list = [x for x in signature_gene_list if x in signature_dict.keys()]

        # Add root marker genes to marker_dict
        for gene_no in range(len(marker_gene_list)):
            marker_dict[marker_gene_list[gene_no]]['Outgroup'] = [root_gene_seqs[gene_no], None, None, None, None, None]
        
        # Get total length of all genes for cluster
        cluster_len_best_list = [int(gene_len_dict.get(x)) for x in marker_gene_list + signature_gene_list]
        cluster_len_best = sum(cluster_len_best_list)

        # Get length of all non-N for each sample 
        cluster_len = [{k: len(v[0].replace('N','')) for k, v in pileup_dict_gene.get(x).items()} for x in marker_gene_list + signature_gene_list if pileup_dict_gene.get(x) is not None]
        cluster_len_samp = {x: sum(filter(None, [y.get(x) for y in cluster_len])) for x in sample_names}
       
        # Select on number of signature and marker genes with at least one read (present in pileup_dict_gene)
        marker_count = [[k for k,v in pileup_dict_gene.get(x).items() if len(v[0])>0] for x in marker_gene_list if pileup_dict_gene.get(x) is not None]
        marker_count_samp = {x: sum([1 for y in marker_count if x in y]) for x in sample_names}
        marker_count_too_few = set([x for x in marker_count_samp if marker_count_samp[x] < min_marker_genes])        

        signature_count = [[k for k,v in pileup_dict_gene.get(x).items() if len(v[0])>0] for x in signature_gene_list if pileup_dict_gene.get(x) is not None]
        signature_count_samp = {x: sum([1 for y in signature_count if x in y]) for x in sample_names}
        signature_count_too_few = set([x for x in signature_count_samp if signature_count_samp[x] < min_signature_genes])        

        # Remove if too many N's
        samp_remove_N = set([k for k,v in cluster_len_samp.items() if v/cluster_len_best < min_nonN])
        samp_remove_count = marker_count_too_few.union(signature_count_too_few)
        
        # Collect samples to remove
        samp_remove = samp_remove_N.union(samp_remove_count)
       
        # Coverage function
        # Sum of depth for all bases divided by sum of sequencing length
        def cov_fun(k, ll):
            ll = [ll.get(x) for x in ll if ll.get(x).get(k) is not None]
            try:
                return sum(filter(None, [y.get(k)[0] for y in ll]))/sum(filter(None, [y.get(k)[1] for y in ll]))
            except ZeroDivisionError:
                return 0

        # Get coverage data for signature genes
        signature_gene_coverage = {x: {k: (v[3],v[5]) for k, v in pileup_dict_gene.get(x).items()} for x in signature_gene_list if pileup_dict_gene.get(x) is not None}
        signature_samp_coverage = {x: cov_fun(x, signature_gene_coverage) for x in sample_names}

        # Get coverage data for marker genes
        marker_gene_coverage = {x: {k: (v[3],v[5]) for k, v in pileup_dict_gene.get(x).items()} for x in marker_gene_list if pileup_dict_gene.get(x) is not None}
        marker_samp_coverage = {x: cov_fun(x, marker_gene_coverage) for x in sample_names}
        
        # Get allele frequency data for signature genes
        signature_gene_af = {x: {k: v[4] for k, v in pileup_dict_gene.get(x).items()} for x in signature_gene_list if pileup_dict_gene.get(x) is not None}
        signature_samp_af = {x: sum(filter(None, [signature_gene_af.get(y).get(x) for y in signature_gene_af])) for x in sample_names}

        # Get allele frequency data for marker genes
        marker_gene_af = {x: {k: v[4] for k, v in pileup_dict_gene.get(x).items()} for x in marker_gene_list if pileup_dict_gene.get(x) is not None}
        marker_samp_af = {x: sum(filter(None, [marker_gene_af.get(y).get(x) for y in marker_gene_af])) for x in sample_names}
        
        # Print stats per cluster
        with open(stat_out, 'a') as fh:
            for i in sample_names:
                fh.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(str(clust),
                                str(i),
                                str(cluster_len_samp[i]/cluster_len_best),
                                str(signature_count_samp[i]),
                                str(signature_samp_coverage[i]),
                                str(signature_samp_af[i]),
                                str(marker_count_samp[i]),
                                str(marker_samp_coverage[i]),
                                str(marker_samp_af[i])))

        # Print stats per gene
        cluster_len_dict = {x: {k: len(v[0].replace('N','')) for k, v in pileup_dict_gene.get(x).items()} for x in signature_gene_list + marker_gene_list if pileup_dict_gene.get(x) is not None}
        nonn_gene = {x: {y: (int(cluster_len_dict.get(y).get(x))/int(gene_len_dict.get(y)) if cluster_len_dict.get(y).get(x) is not None else 0) for y in cluster_len_dict} for x in sample_names}
        
        with open(stat_gene_out, 'a') as fh:
            # Per sample
            for i in sample_names:
                # Per gene
                gene_sub = nonn_gene[i]
                for j in gene_sub:
                    
                    # Get multiallelic fraction
                    AF_sig = af_gene(signature_gene_af, gene_len_dict, i, j)
                    AF_mark = af_gene(marker_gene_af, gene_len_dict, i, j)
                    
                    # Get coverage info
                    sig_cov = cov_gene(signature_gene_coverage, i, j)
                    mark_cov = cov_gene(marker_gene_coverage, i, j)

                    fh.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        str(clust),
                        str(i),
                        str(j),
                        str(gene_sub[j]),
                        str(sig_cov),
                        str(mark_cov),
                        str(AF_sig),
                        str(AF_mark)))


        # Delete samples
        marker_dict = delete_keys_dict(marker_dict, samp_remove)
        signature_dict = delete_keys_dict(signature_dict, samp_remove)
        
        # Remove previous alignments
        try:
            os.mkdir(os.path.join(align_out, clust))
        except FileExistsError:
            shutil.rmtree(os.path.join(align_out, clust))
            os.mkdir(os.path.join(align_out, clust))
        
        # Write to fastas
        for gene_id in marker_dict:
            sub_dct = marker_dict[gene_id]
            with open(os.path.join(align_out, clust, clust+'_'+gene_id+'.fna'), 'w') as fh:
                for fa in sub_dct:
                    if len(sub_dct[fa][0]) > 0:
                        fh.write('>'+fa+'\n'+sub_dct[fa][0]+'\n')

        for gene_id in signature_dict:
            sub_dct = signature_dict[gene_id]
            with open(os.path.join(align_out, clust, clust+'_'+gene_id+'.fna'), 'w') as fh:
                for fa in sub_dct:
                    if len(sub_dct[fa][0]) > 0:
                        fh.write('>'+fa+'\n'+sub_dct[fa][0]+'\n')

def af_gene(dd, gdd, samp, gene):
    '''
    With a nested dictionary (dd) of multiallelic counts for genes (first layer) and samples (second layer),
    and a dictionary (gdd) of gene lengths
    return the multiallelic site fraction for the gene (gene) for that sample (samp)
    '''
    try:
        AF = dd.get(gene).get(samp)
        try:
            AF = AF / int(gdd.get(gene))
        except:
            AF = 0
    except:
        AF = 'NA'
    return(AF)

def cov_gene(dd, samp, gene):
    '''
    With a nested dictionary (dd) of total depth (3) and seq_len (5) for genes (first layer) and samples (second layer) 
    return the coverage for the gene (gene) for that sample (samp)
    '''
    try:
        cov = dd.get(gene).get(samp)
        if cov is None:
            cov = 0
        else:
            try:
                cov = cov[0]/cov[1]
            except ZeroDivisionError:
                cov = 0
    except:
        cov = 'NA'
    return(cov)

def delete_keys_dict(dct, keys):
    '''
    Will delete keys (keys) from the second layer dictionaries of a deep copy of the input dictionary (dct)
    '''
    dct_copy = copy.deepcopy(dct)
    for k in keys:
        for kk in dct:
            try:
                del dct_copy[kk][k]
            except KeyError:
                pass
    return(dct_copy)

def inverse_dict(dd):
    '''
    Reverse the inner and outer dicts of a nested dictionary
    such that dd[samp][gene] becomes dd[gene][samp]
    '''
    inversed_dict = defaultdict(dict)
    for k, v in dd.items():
        for sk, sv in v.items():
            inversed_dict[sk][k] = sv

    return(dict(inversed_dict))


def parse_base_pairs(ref_base, x, min_af, min_depth):
    '''
    Parse mpileup read base column
    
    ref_base: Base of reference sequences
    x: pileup read base string
    min_af: Minimum allele frequency, if below this an N is returned

    Return:
    Str: Base
    Int: Multiallelic sites
    '''
   
    # Allele frequency indicator
    af_low = 0

    # Ignore positions of reads
    x = re.sub('\^.', '', x)
    x = x.replace('$', '')

    # Force reverse reads to forward
    x = x.upper()
    x = x.replace(',', '.')

    # Insert reference
    x = x.replace('.', ref_base)
    
    # Return N if depth too low
    if len(x) < min_depth:
        seq = 'N'
        return seq, af_low

    # If no indels, just report most seen and the frequency
    if (x.find('+') < 0) and (x.find('-') < 0):
        chars = ('*','A','C','T','G','N')
        counts = [x.count(y) for y in chars]
        seq = chars[counts.index(max(counts))]
        
        # If deletion is most common, insert no sequence
        if seq == '*':
            seq = ''
   
        # If allele frequency is low, replace with N
        if max(counts)/len(x) < min_af:
            seq = 'N'
            af_low = 1

    # If indels present, traverse through string
    else:

        # Initialize variables
        char_dict = {}
        deletion = False
        insertion = False
        del_count = 0
        ins_count = 0
        ins_finished = False
        s_prev = ''

        for s in x:
            # If in a deletion, skip this char
            if del_count > 0:
                del_count -= 1
                continue

            # If in an insertion, combine characters
            if ins_count > 0:
                ins_count -= 1
                ins_seq = ins_seq + s
                
                # When insertion finished, overwrite the s var
                if ins_count == 0:
                    s = s_prev + ins_seq
                    ins_finished = True
                else:
                    continue

            # If minus was detected in previous char, grap the quantifier
            if deletion:
                del_count = int(s)
                deletion = False
                continue

            # If plus was detected in previous char, grap the quantifier
            if insertion:
                ins_count = int(s)
                insertion = False
                ins_seq = ''
                continue 
            
            # If indel, enter new states
            if s == '-':
                deletion = True
            elif s == '+':
                insertion = True

            # If not in indel save in dict
            else:
                if s in char_dict:
                    char_dict[s] += 1
                else:
                    char_dict[s] = 1

                # If end of insertion, remove count of previous
                if ins_finished:
                    char_dict[s_prev] -= 1
                    ins_finished = False
                
                # Save previous char for insertions
                s_prev = s

        counts = list(char_dict.values())
        seq = list(char_dict.keys())[counts.index(max(counts))]
        
        # If deletion most common, insert no sequence
        if seq == '*':
            seq = ''

        # If allele frequency low, replace with N
        if max(counts)/len(x) < min_af:
            seq = 'N'
            af_low = 1

    return seq, af_low


def mpileup_read(file, min_af = const_min_af, min_depth = const_min_depth):
    '''
    Read a samtools mpileup file and output sequence dict

    min_af: Minimum allowed allele frequency, if above this an N is returned at that position

    Return:
    dict: key, value = sequence ID, [sequence, first position, last position, total depth, n bases with too low allele freq, seq_len]
    
    '''

    with open(file, 'r') as file_handle:

        # Initialize dict
        seq_dict = {}

        # Initialize variables for first line
        seq_id = None
        first = 0
        last = 0
        depth = 0
        af_count = 0
        seq_len = 0

        # Read data
        for line in file_handle:
                
            line_data = line.rstrip().split('\t')
            
            # Starting with new reference
            if line_data[0] != seq_id:
               
                # If this is first line in file, don't save
                if seq_id is None:
                    pass
                
                # If not first line in file, save results
                else:
                    seq_dict[seq_id] = [sequence, first, last, depth, af_count, seq_len]
                
                # Start variables for new reference
                seq_id = line_data[0]
                first = 0
                last = 0
                
                if int(line_data[3]) == 0:
                    seq_base = ''
                    af_ind = 0
                else:
                    seq_base, af_ind = parse_base_pairs(line_data[2], line_data[4], min_af, min_depth)
                
                sequence = seq_base
                af_count = af_ind

                seq_len = len(seq_base)
                depth = int(line_data[3])*len(seq_base)
         
                # Only save first if not a deletion
                if len(seq_base) > 0: 
                    first = int(line_data[1])
                else:
                    first = 0

            # Continuing with same reference as previous line
            else:
                if int(line_data[3]) == 0:
                    seq_base = ''
                    af_ind = 0
                else: 
                    seq_base, af_ind = parse_base_pairs(line_data[2], line_data[4], min_af, min_depth)
                
                sequence += seq_base
                af_count += af_ind

                seq_len += len(seq_base)
                depth += int(line_data[3])*len(seq_base)

                # If this is the first character
                if first == 0:
                    if len(seq_base) > 0:
                        first = int(line_data[1])

                # Save base if not deletion
                if len(seq_base) > 0:
                    last = int(line_data[1])

        # End of file, save results
        if len(sequence) > 0: 
            seq_dict[seq_id] = [sequence, first, last, depth, af_count, seq_len]

        return seq_dict



def read_refs(ref_fasta):
    '''
    Load the reference sequences to get the lengths
    '''
    len_dict = {}
    with open(ref_fasta, 'r') as fh:
        for fa in SeqIO.parse(fh, 'fasta'):
            len_dict[fa.id] = len(fa.seq)

    return(len_dict)

def read_roots(root_fasta):
    '''
    Load the root marker gene sequences
    '''
    seq_dict = {}
    with open(root_fasta, 'r') as fh:
        for fa in SeqIO.parse(fh, 'fasta'):
            seq_dict[fa.id] = str(fa.seq)

    return(seq_dict)

#######################################################################################################################
if __name__ == '__main__':
    main()
