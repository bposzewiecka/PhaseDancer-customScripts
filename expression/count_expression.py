import csv
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from Bio import SeqIO
import pysam

def get_sample_size(fn):

    if fn is None: return 0

    sample_size = 0
    
    for record in SeqIO.parse(fn, "fastq"):
        sample_size += len(record.seq)

    return sample_size

def count_expression(data_run, transcripts_fn, bam_fn, fasta1_fn, fasta2_fn, output_fn):
    
    sample_size = get_sample_size(fasta1_fn) + get_sample_size(fasta2_fn)

    exons_by_gene = defaultdict(IntervalTree)
    genes = {}

    with open(transcripts_fn) as csvfile:

        reader = csv.DictReader(csvfile,  delimiter='\t')

        for row in reader:

            start = int(row['chromStart'])
            end = int(row['chromEnd'])

            for block_size, chrom_start in zip(map(int, row['blockSizes'].split(',')[:-1]), map(int, row['chromStarts'].split(',')[:-1])):   

                exon =  Interval(start + chrom_start,  start + chrom_start + block_size)
                exons_by_gene[row['geneName']].add(exon)

            if row['geneName'] in genes:

                if genes[row['geneName']]['start'] > start:
                    genes[row['geneName']]['start'] =  start
                if genes[row['geneName']]['end'] <  end:
                    genes[row['geneName']]['end'] = end

            else:
                genes[row['geneName']] = { 'chrom': row['#chrom'], 'start':  start, 'end': end }

    samfile = pysam.AlignmentFile(bam_fn, "rb")

    with open(output_fn, "w") as f_out:    
        
        for gene in genes:
            exons_by_gene[gene].merge_overlaps()    

        for gene, exons  in exons_by_gene.items(): 

            expression = 0
            exons_size = 0

            for exon in exons:

                for pileupcolumn in samfile.pileup(genes[gene]['chrom'], exon.begin, exon.end, truncate=True):        
                    expression += pileupcolumn.n

                exons_size += exon.end - exon.begin

            genes[gene]['expression'] = expression 
            genes[gene]['exons_size'] = exons_size
            genes[gene]['sample_size'] = sample_size
            genes[gene]['normalised_expression'] = expression / exons_size / sample_size * 1000 * 1000 * 1000

            f_out.write('\t'.join(map(str, data_run + [ gene, expression, exons_size, sample_size, genes[gene]['normalised_expression']])))
            f_out.write('\n')

    samfile.close()
