import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def load_data(antibody, evo_round):
    fname = f'data/rounds/{antibody}_{evo_round}.txt'
    data = []
    with open(fname) as f:
        f.readline()
        for line in f:
            fields = line.rstrip().split('\t')
            antigen = fields[0]
            chain = fields[1]
            mut_str = fields[2]
            try:
                fold_change = float(fields[-1])
                log_fold_change = np.log10(float(fields[-1]))
            except:
                fold_change = 0.
                log_fold_change = float('nan')
            data.append([
                antibody,
                evo_round,
                antigen,
                chain,
                mut_str,
                fold_change,
                log_fold_change,
            ])
    return data

            
if __name__ == '__main__':
    
    antibodies = [
        'medi',
        'uca',
        'mab114',
        'mU',
        's309',
        'r7',
        'c143',
    ]

    rounds = [
        'alpha',
        'beta'
    ]

    data = []
    for antibody in antibodies:
        for evo_round in rounds:
            data += load_data(antibody, evo_round)

    df = pd.DataFrame(data, columns=[
        'antibody',
        'evo_round',
        'antigen',
        'chain',
        'mut_str',
        'fold_change',
        'log_fold_change',
    ])

    plt.figure()
    sns.catplot(
        data=df,
        x='evo_round',
        y='fold_change',
        row='antibody',
        sharey=False,
        s=20,
        jitter=0.15,
        aspect=1.1,
    )
    plt.savefig('figures/rounds.svg')
    plt.close()

    plt.figure(figsize=(5, 20))
    sns.catplot(
        data=df,
        x='evo_round',
        y='log_fold_change',
        row='antibody',
        sharey=False,
    )
    plt.savefig('figures/rounds_log_scale.svg')
    plt.close()
    
