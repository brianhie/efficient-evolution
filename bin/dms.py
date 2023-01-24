import numpy as np
import pandas as pd
import scipy.stats as ss

from amis import reconstruct_multi_models

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='DMS sequence analysis')
    parser.add_argument(
        '--model-names',
        type=str,
        default=[ 'esm1b', 'esm1v1', 'esm1v2', 'esm1v3', 'esm1v4', 'esm1v5', ],
        nargs='+',
        help='Type of language model (e.g., esm1b, esm1v1)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=None,
        help='Softness of reconstruction'
    )
    args = parser.parse_args()
    return args

def dms_results(fname, model_names, alpha=None,):
    namespace = fname.split('/')[-1].split('.')[0]
    with open(fname):
        df = pd.read_csv(fname, delimiter=',')

    wt_seq = ''
    for mut in df['variant']:
        aa_orig, pos = mut[0], int(mut[1:-1]) - 1
        if len(wt_seq) > pos:
            assert(wt_seq[pos] == aa_orig)
        else:
            wt_seq += aa_orig

    mutations_models = reconstruct_multi_models(wt_seq, model_names, alpha=alpha)

    mutation_scores, studies = {}, []
    for column in df.columns:
        if column.startswith('DMS'):
            studies.append(column)
            for mut, val in zip(df['variant'], df[column]):
                pos, aa_mut, aa_wt = int(mut[1:-1]) - 1, mut[-1], mut[0]
                assert(wt_seq[pos] == aa_wt)
                key = (pos, aa_wt, aa_mut)
                if key not in mutation_scores:
                    mutation_scores[key] = {}
                mutation_scores[key][column] = val

    for study in studies:
        for mutation, n_models in sorted(
                mutations_models.items(), key=lambda item: item[1]
        ):
            assert(mutation in mutation_scores)
            fields = [
                mutation,
                mutation_scores[mutation][study],
                n_models,
                namespace,
                study,
            ]
            print('\t'.join([ str(field) for field in fields ]))

if __name__ == '__main__':
    args = parse_args()

    fnames = [
        'data/dms/dms_adrb2.csv',
        'data/dms/dms_bla.csv',
        'data/dms/dms_env.csv',
        'data/dms/dms_ha_h1.csv',
        'data/dms/dms_ha_h3.csv',
        'data/dms/dms_infa.csv',
        'data/dms/dms_mapk1.csv',
        'data/dms/dms_p53.csv',
        'data/dms/dms_pafa.csv',
    ]

    for fname in fnames:
        dms_results(fname, args.model_names, alpha=args.alpha)
