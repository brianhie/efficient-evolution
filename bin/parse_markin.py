from Bio import SeqIO
import pandas as pd
import sys

from utils import deep_mutational_scan

if __name__ == '__main__':
    df = pd.read_csv('data/pafa/abf8761_markin_data-s1.csv')
    wt_seq = str(SeqIO.read('data/pafa/pafa_wt_seq.fa', 'fasta').seq)
    offset = 20

    mut_pheno = {}
    for variant, kcat, km, kcat_km in zip(
            df['variant'],
            df['kcat_cMUP_s-1'],
            df['KM_cMUP_uM'],
            df['kcatOverKM_cMUP_M-1s-1'],
    ):
        if '/' in variant or variant == 'WT':
            continue
        wt, mt, pos = variant[0], variant[-1], int(variant[1:-1]) - offset
        mut_str = f'{wt}{pos + 1}{mt}'
        assert(wt_seq[pos] == wt)
        mut_pheno[mut_str] = {
            'kcat': kcat,
            'km': km,
            'kcat_km': kcat_km,
        }

    header = [
        'variant',
        'DMS_kcat',
        'DMS_km',
        'DMS_kcat_km',
    ]
    print(','.join(header))
        
    for pos, wt, mt in deep_mutational_scan(wt_seq):
        mut_str = f'{wt}{pos + 1}{mt}'
        fields = [ mut_str ]
        for field in header[1:]:
            field = field[4:]
            if mut_str in mut_pheno:
                fields.append(mut_pheno[mut_str][field])
            else:
                fields.append('')
        print(','.join([ str(field) for field in fields ]))
