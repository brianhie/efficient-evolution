from Bio import SeqIO
import sys

seq_fasta = open('data/uniref/uniref202003/uniref_2020_03.fasta', 'w')
seq_mapping = open('data/uniref/uniref202003/uniref_2020_03_mapping.txt', 'w')

in_sequence, in_member = False, False
n_seqs = 0

with open('data/uniref/uniref202003/uniref90.xml') as f:
    for line in f:
        if line.startswith('<sequence '):

            if line.rstrip().endswith('</sequence>'):
                curr_seq = line.split('>')[1].split('<')[0]
                n_seqs += 1
                seq_fasta.write(f'>seq{n_seqs}\n')
                seq_fasta.write(f'{curr_seq}\n')
                continue
            
            in_sequence = True
            curr_seq = ''
            n_seqs += 1
            continue

        if line.startswith('</sequence>'):
            in_sequence = False
            seq_fasta.write(f'>seq{n_seqs}\n')
            seq_fasta.write(f'{curr_seq}\n')
            continue

        if line.startswith('<member>'):
            in_member = True

        if line.startswith('</member>'):
            in_member = False

        if in_sequence:
            curr_seq += line

        if in_member:
            if 'UniProtKB ID' in line or \
               'UniProtKB accession' in line or\
               'UniParc ID' in line:
                value = line.rstrip().rstrip('"/>').split('"')[-1]
                seq_mapping.write('\t'.join([
                    value,
                    f'seq{n_seqs}'
                ]) + '\n')

seq_fasta.close()
seq_mapping.close()
