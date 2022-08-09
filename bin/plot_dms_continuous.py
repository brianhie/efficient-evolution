import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as ss
import seaborn as sns
from sklearn.metrics import average_precision_score as aupr

np.random.seed(0)


def load_dms_lm(fname):
    data = []
    with open(fname) as f:
        for line in f:
            fields = line.rstrip().split()
            if fields[3] == 'nan':
                continue
            mut_str = ' '.join(fields[:3])
            fitness = float(fields[3])
            n_models = int(fields[4])
            namespace = fields[5]
            dms_name = fields[6]
            data.append([
                namespace,
                dms_name,
                mut_str,
                n_models,
                fitness,
            ])

    return pd.DataFrame(data, columns=[
        'namespace',
        'dms_name',
        'mut_str',
        'n_models',
        'fitness',
    ])


def pct_norm(array, fitness_min, fitness_max):
    return (array - fitness_min) / (fitness_max - fitness_min)


if __name__ == '__main__':
    settings = [
        [ 'dms_adrb2', 'DMS_0.625', 2.8, 9 ],
        [ 'dms_bla', 'DMS_amp_2500_(b)', 0.01, 10 ],
        [ 'dms_env', 'DMS', 0.1, 31 ],
        [ 'dms_ha_h1', 'DMS', 0.1, 32 ],
        [ 'dms_ha_h3', 'DMS', 0.1, 16 ],
        [ 'dms_infa', 'DMS_min', 0.98, 10 ],
        [ 'dms_mapk1', 'DMS_SCH', 2.5, 13 ],
        [ 'dms_p53', 'DMS_null_etoposide', 1, 17 ],
    ]

    data, random = [], []
    for namespace, dms_name, cutoff, sample_size in settings:
        df = pd.read_csv(f'data/dms/{namespace}_esm.csv', delimiter=',')
        n_positive = sum(df[dms_name] > cutoff)
        n_total = sum(np.isfinite(df[dms_name]))
        background_frac = n_positive / n_total
        random.append([ namespace, background_frac * sample_size, ])

        for column in df.columns:
            if column == 'variant' or column.startswith('DMS'):
                continue
            if sum(np.isfinite(df[column])) == 0:
                continue
            if 't33_650M_' in column:
                continue

            prediction = df[column].values
            experiment = df[dms_name].values
            n_positives = []
            for i in range(100): # Break ties with random guessing.
                rand_idx = np.random.permutation(len(prediction))
                prediction = prediction[rand_idx]
                experiment = experiment[rand_idx]
                idx_topk = np.argsort(-prediction)[:sample_size]
                n_positives.append(sum((experiment > cutoff)[idx_topk]))
            n_positive_topk = np.mean(n_positives)

            finite_idx = np.logical_and(np.isfinite(prediction), np.isfinite(experiment))
            prediction_finite = prediction[finite_idx]
            experiment_finite = experiment[finite_idx]

            data.append([
                namespace,
                column,
                n_positive_topk,
                # Below not really applicable to ESM_vote but compute anyway.
                abs(ss.spearmanr(experiment_finite, prediction_finite)[0]),
                aupr(experiment_finite > cutoff, prediction_finite),
            ])
            
    df = pd.DataFrame(data, columns=[
        'namespace',
        'model_name',
        'n_positive_topk',
        'spearmanr',
        'aupr',
    ])

    df.to_csv('dms_continuous.tsv', sep='\t')

    for metric in [ 'n_positive_topk', 'spearmanr', 'aupr', ]:
        namespaces = sorted(set(df['namespace']))
        model_ranks = {}
        for name in namespaces:
            df_namespace = df[df['namespace'] == name]
            for model, rank in zip(
                    df_namespace['model_name'],
                    ss.rankdata(-df_namespace[metric], method='average')
            ):
                model = model.split(' ')[0]
                if model not in model_ranks:
                    model_ranks[model] = []
                model_ranks[model].append(rank)
        df_mean_rank = pd.DataFrame([
            [ model, np.mean(model_ranks[model]), len(model_ranks[model]) ]
            for model in model_ranks
        ], columns=[ 'model', 'mean_rank', 'n_ranks', ])
        df_mean_rank.sort_values('mean_rank', inplace=True)
        df_mean_rank.to_csv(f'dms_continuous_mean_rank_{metric}.tsv', sep='\t')
    
        top_models = df_mean_rank.head(6)
        to_highlight = {
            model for idx, model in enumerate(top_models['model'])
        }
        colors = [
            '#e8bcf9',
            '#fdda0d',
            '#add8e6',
            '#afe1af',
            '#ffa500',
            '#c41e3a',
        ]
        plt.figure()
        sns.stripplot(
            data=df[[ name not in to_highlight for name in df['model_name'] ]],
            x='namespace',
            y=metric,
            color='#aaaaaa',
            size=6,
        )
        sns.swarmplot(
            data=df[[ name in to_highlight for name in df['model_name'] ]],
            x='namespace',
            y=metric,
            hue='model_name',
            palette=colors,
            size=7,
        )
        if metric == 'n_positive_topk':
            sns.boxplot( # Plot random guessing as dashed line.
                showmeans=True,
                meanline=True,
                meanprops={'color': 'k', 'ls': '--', 'lw': 1},
                medianprops={'visible': False},
                whiskerprops={'visible': False},
                zorder=10,
                x='namespace',
                y='n_positive_topk',
                data=pd.DataFrame(random, columns=[ 'namespace', 'n_positive_topk' ]),
                showfliers=False,
                showbox=False,
                showcaps=False,
            )
        elif metric == 'spearmanr':
            plt.axhline(y=0., c='k', linestyle='--',)
        plt.savefig(f'figures/dms_continuous_{metric}.svg')
        plt.close()
