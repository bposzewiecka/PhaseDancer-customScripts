import os
import edlib
import re, math
from Bio import SeqIO
from collections  import Counter

from src.utils.yaml_utils import save_yaml, load_yaml
from src.utils.configuration import snakemake_validate_config, CONFIG_FN
#from src.scripts.clustering import AlignmentArray
#from src.scripts.cluster_selecting import read_clusters

import sys

to_trim = 'CANUBAC-0036-AC275606', 'CANUBAC-0057-AC278733', 'CANUBAC-0060-AC278342', 'CANUBAC-0062-AC278737', 'CANUBAC-0081-AC278480', 'CANUBAC-0091-AC278602', 'CANUBAC-0105-AC278308', 'CANUBAC-0118-AC278544', 'CANUBAC-0119-AC278687', 'CANUBAC-0141-AC278751', 'CANUBAC-0148-AC278524', 'CANUBAC-0183-AC278803', 'CANUBAC-0232-AC278630', 'CANUBAC-0284-AC278412', 'CANUBAC-0294-AC278439', 'CANUBAC-0295-AC278652', 'CANUBAC-0323-AC278861'

def get_sequence(fasta_in_fn, trim = False):
   
    for record in SeqIO.parse(fasta_in_fn, 'fasta'):
        trim_size = 0
        if trim:
            trim_size = 500

        return record.name, str(record.seq)[trim_size:] 

def get_identity_by_cigar(cigar):

    alignment_size = 0
    errors = 0
    all_bases = 0

    for bases_sign, sign in re.findall('(\d+(=|D|I|X))', cigar):
        bases = int(bases_sign[:-1])

        alignment_size += bases

        if sign not in ('=', 'X'):
           errors += bases
        
        all_bases += bases

    print(cigar)

    return {
        'alignment_size': alignment_size,
        'edit_distance': errors,
        'identity': 1 - errors / alignment_size,
        'phred_quality': math.log10(errors/alignment_size) * -10 if errors else None
    }

def get_alignment_results(reference, assembly):

    result = edlib.align(reference,assembly, mode="HW", task="path")

    return get_identity_by_cigar(result['cigar'])
       

def align_assembled_contig(contig_fn, assembled_contig_fn):
   
    contig_name, contig = get_sequence(contig_fn)
    _, assembled_contig = get_sequence(assembled_contig_fn)

    alignment = edlib.align(contig, assembled_contig, mode="HW", task="path")

    return alignment

PHASEDANCER_DATA_DIR = os.environ['PHASEDANCER_DATA_DIR'] 
OUTPUT_DIR = PHASEDANCER_DATA_DIR + '/'

configfile: OUTPUT_DIR + CONFIG_FN

errors = snakemake_validate_config(OUTPUT_DIR + CONFIG_FN, check_files = False)

def get_config_value(key, sample, contig = None):

    value = None

    if key in config['samples'][sample]:
        value = config['samples'][sample][key]
	
    if contig is not None and type(config['samples'][sample]['contigs']) == dict and key in config['samples'][sample]['contigs'][contig]:
        value = config['samples'][sample]['contigs'][contig][key]

    return value

def get_params(sample,  contig):

    return { 'assembler': 'minimap2', 'cluster': '000', 'start_number': '000', 'sample': sample, 'contig': contig, 'end_number':  get_config_value('iterations', sample, contig), 'cl_type': 'nc'}

def get_input_dicts():

    def get_sample_params(sample):
        
        contigs = get_config_value('contigs', sample)
	
        if type(contigs) == dict:

            contigs = contigs.keys()

        return [ get_params(sample, contig) for contig in contigs]

    if config['sample'] == 'all':

        data = []
  
        for sample in config['samples']:
            data += get_sample_params(sample)
	
        return data

    elif config['contig'] == 'all':

         return get_sample_params(config['sample'])

    else:

        return [ get_params(config['sample'],  config['contig'])]
	
def is_assembly_successfull(fn):
    try:
        with open(fn) as f:
           l = f.readline()
           return l[0] == '0'
    except:
        return False

rule main:
    input:
        [ expand(OUTPUT_DIR + 'data/{sample}/{contig}/{contig}.status.txt', **input_dict)  for input_dict  in get_input_dicts()],

rule run_assembly: 
    input:
        OUTPUT_DIR + 'data/{sample}/{contig}/{contig}.fasta'
    output:
        status = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}.status.txt'
    params:
        assembly_log = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}.log.txt'
    threads:
        3
    shell:
        """
        set +e
        snakemake --config sample={wildcards.sample} contig={wildcards.contig} --cores {threads} --rerun-incomplete >> {params.assembly_log}
        exitcode=$?

        echo $exitcode > {output.status}
        exit 0
        """

rule run_statistics:
    input:
        status = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}.status.txt',
        assembly = lambda wildcards: expand(OUTPUT_DIR + 'data/{sample}/{contig}/{assembler}/{cl_type}_{cluster}/merged_contigs/seq_{start_number}.seq_{end_number}.{contig}.{cl_type}_{cluster}.{sample}.{assembler}.merged_contig.fasta', **get_params(wildcards.sample, wildcards.contig))[0],
        reference =  OUTPUT_DIR + 'data/{sample}/{contig}/{contig}_bac.fasta'
    output:
        stats = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}.assembly_stats.yaml'
    run:
        trim = wildcards.contig in to_trim
        _, reference = get_sequence(input.reference, trim = trim)
        _, assembly =  get_sequence(input.assembly)    
                  
        results = get_alignment_results(reference, assembly)
	print(results)
        save_yaml(output.stats, results)


rule run_statistics2:
    input:
        bac = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}_bac.fasta'
    output:
        stats = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}.contigs_on_bac.yaml'
    run:
        bac_fn = input.bac

        sample = wildcards.sample
        contig = wildcards.contig

        results = { 'results' : []}

        try:
            for no in range(1000):

                contig_fn = OUTPUT_DIR + f'data/{sample}/{contig}/minimap2/nc_000/seq_{no:03d}/seq_{no:03d}.{contig}.nc_000.{sample}.minimap2.uncompressed.fasta'
                             
                alignment =  align_assembled_contig(contig_fn, bac_fn)
            
                results['results'].append( {"no": no, "result":  alignment })

        except Exception as e:
            print(e)

        save_yaml(output.stats, results)


rule run_statistics3:
    input:
        bac = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}_bac.fasta'
    output:
        stats = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}.contigs_on_assembled_contigs.yaml'
    run:
        sample = wildcards.sample
        contig = wildcards.contig

        results = { 'results' : []}

        try:
            for no in range(1000):

                contig_fn = OUTPUT_DIR + f'data/{sample}/{contig}/minimap2/nc_000/seq_{no:03d}/seq_{no:03d}.{contig}.nc_000.{sample}.minimap2.uncompressed.fasta'
                assembled_contig_fn = OUTPUT_DIR + f'data/{sample}/{contig}/minimap2/nc_000/seq_{no:03d}/seq_{no:03d}.{contig}.nc_000.{sample}.minimap2.assembled_contig.fasta'
 
                alignment = align_assembled_contig(contig_fn, assembled_contig_fn)

                results['results'].append( {"no": no, "result":  alignment })

        except:
            pass

        save_yaml(output.stats, results)

rule cluster_errors:
    input:
        bac = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}_bac.fasta'
    output:
        stats = OUTPUT_DIR + 'data/{sample}/{contig}/{contig}.cluster_errors.yaml'
    run: 

        sample = wildcards.sample
        contig = wildcards.contig

        results = { 'results' : []}

        prev_errors = 0
        curr_errors = 0

        try:
            for no in range(1,100):

                selected_clusters_fn = OUTPUT_DIR + f'data/{sample}/{contig}/minimap2/nc_000/seq_{no:03d}/clusters/seq_{no:03d}.{contig}.nc_000.{sample}.minimap2.selected_cluster.yaml'
		clusters_fn =  OUTPUT_DIR + f'data/{sample}/{contig}/minimap2/nc_000/seq_{no:03d}/clusters/seq_{no:03d}.{contig}.nc_000.{sample}.minimap2.clusters.tsv'

		bam_fn = OUTPUT_DIR + f'data/{sample}/{contig}/minimap2/nc_000/seq_{no:03d}/seq_{no:03d}.{contig}.nc_000.{sample}.minimap2.compressed.bam'
                fasta_fn = OUTPUT_DIR + f'data/{sample}/{contig}/minimap2/nc_000/seq_{no:03d}/seq_{no:03d}.{contig}.nc_000.{sample}.minimap2.compressed.fasta'

                clusters_info = load_yaml(selected_clusters_fn)

                prev_selected = clusters_info['prev_selected'] 
                selected = clusters_info['selected'] 

		cluster = read_clusters(clusters_fn)[selected]

                cluster_stats = clusters_info['cluster_stats']
                prev_stats =  [v for v in cluster_stats[prev_selected].values() if v > 2 ] 
                prev_errors += len(prev_stats) > 1

                curr_stats = [ v[selected]  for k, v in cluster_stats.items() if v[selected]  ] 
                curr_errors += len(curr_stats) > 1
		



                if len(curr_stats) > 1:
                    alignment_array = AlignmentArray(bam_fn, fasta_fn)

		    insertions_by_ref = alignment_array.insertions_by_ref_pos

                    diff = False 

                    for ref_pos  in sorted(insertions_by_ref.keys()):
 
                        insertions = insertions_by_ref[ref_pos]
			#print(ref_pos)
			#print(insertions)
                        in_cluster = [ v  for k, v in insertions.items() if k in cluster ] 
			#print(in_cluster)
		
                        c = Counter(in_cluster)
      
                        if len(c) > 0 and len(in_cluster)>2:
                            print(ref_pos, c)
                            diff = True

                    if diff: print('Diff')
		 	
                    print('Cluster:', no, wildcards.contig, len(cluster))
		
	
                results['results'].append({
                    'seq': no,
                    'prev_selected':  prev_selected, 
                    'selected':  selected, 
                    'prev_stats': prev_stats,
                    'curr_stats': curr_stats,
                })

        except Exception as e:
            print(e)
            pass
        
        results['curr_errors'] = curr_errors
        results['prev_errors'] = prev_errors
        results['all_errors'] = prev_errors + curr_errors

        save_yaml(output.stats, results)
