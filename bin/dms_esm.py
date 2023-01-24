import numpy as np
import pandas as pd
import scipy.stats as ss

from amis import reconstruct_multi_models
from predict_esm import predict_esm

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='DMS sequence analysis')
    parser.add_argument(
        '--model-locations',
        type=str,
        default=[
            'esm1b_t33_650M_UR50S',
            'esm1v_t33_650M_UR90S_1',
            'esm1v_t33_650M_UR90S_2',
            'esm1v_t33_650M_UR90S_3',
            'esm1v_t33_650M_UR90S_4',
            'esm1v_t33_650M_UR90S_5',
        ],
        nargs='+',
        help='Type of language model (e.g., esm1b, esm1v1)'
    )
    args = parser.parse_args()
    return args

def dms_esm(fname, model_locations):
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

    df_esm = predict_esm(wt_seq, model_locations,
                         mutation_col='variant', offset_idx=1)
    summed_log_probs = []
    for idx, row in df_esm.iterrows():
        summed_log_probs.append(sum([
            row[model_location] for model_location in model_locations
        ]))
    df_esm['ESM_summed'] = summed_log_probs

    model_names = [ 'esm1b', 'esm1v1', 'esm1v2', 'esm1v3', 'esm1v4', 'esm1v5', ] 
    if '_infa' in fname:
        alpha = 0.5
    else:
        alpha = None
    mutations_models, mutations_model_names = reconstruct_multi_models(
        wt_seq,
        models_names,
        alpha=alpha,
        return_names=True,
    )
    esm_vote, esm_vote_component = [], {
        name: [] for _, names in mutations_model_names.items()
        for name in names
    }
    for mutation in df_esm['variant']:
        pos, wt, mt = int(mutation[1:-1]) - 1, mutation[0], mutation[-1]
        esm_vote.append(mutations_models.get((pos, wt, mt), 0))
        model_names = mutations_model_names.get((pos, wt, mt), [])
        for name in esm_vote_component:
            esm_vote_component[name].append(int(name in model_names))
    df_esm['ESM_vote'] = esm_vote
    for name, values in esm_vote_component.items():
        df_esm[name + '_vote'] = values
    
    df.set_index('variant', inplace=True)
    df_esm.set_index('variant', inplace=True)
    df = df.join(df_esm, how='left')
    df.to_csv('.'.join(fname.split('.')[:-1]) + '_esm.csv')

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
        print(fname)
        dms_esm(fname, args.model_locations)
