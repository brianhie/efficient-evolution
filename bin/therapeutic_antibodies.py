import pandas as pd
import time

from amis import reconstruct_multi_models
from utils import tprint


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description='Therapeutic antibody sequence analysis'
    )
    parser.add_argument(
        'seq_fname',
        type=str,
        help='Path to therapeutic antibody file from Thera-SAbDab'
    )
    parser.add_argument(
        'ofname',
        type=str,
        help='File to save results'
    )
    parser.add_argument(
        '--model-names',
        type=str,
        default=[ 'esm1b', 'esm1v1', 'esm1v2', 'esm1v3', 'esm1v4', 'esm1v5', ],
        nargs='+',
        help='Type of language model (e.g., esm1b, esm1v1)'
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    df = pd.read_csv(args.seq_fname, sep=',')

    data = []
    time_start = time.time()
    for idx, (name, vh, vl, vh_bi, vl_bi) in enumerate(zip(
            df['Therapeutic'],
            df['Heavy Sequence'],
            df['Light Sequence'],
            df['Heavy Sequence (if bispec)'],
            df['Light Sequence (if bispec)'],
    )):
        tprint(f'Evaluating {name}')
        if idx % 10 == 0:
            now = time.time()
            tprint(f'{now - time_start} seconds until iteration {idx}')

        for wt_seq in [ vh, vl, vh_bi, vl_bi ]:
            if wt_seq == 'na':
                continue
            mutations_models = reconstruct_multi_models(wt_seq, args.model_names)
            for mutation, n_models in sorted(
                    mutations_models.items(), key=lambda item: item[1]
            ):
                mutation_str = f'{mutation[1]}{mutation[0] + 1}{mutation[2]}'
                data.append([
                    name,
                    wt_seq,
                    mutation_str,
                    n_models,
                ])

    df_results = pd.DataFrame(data, columns=[
        'Therapeutic',
        'Variable Sequence',
        'Mutation',
        'Number of Language Models',
    ])
    df_results.to_csv(args.ofname)
