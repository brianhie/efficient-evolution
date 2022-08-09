import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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
        [ 'dms_adrb2', 'DMS_0.625', 2.8, 2 ],
        [ 'dms_bla', 'DMS_amp_2500_(b)', 0.01, 2 ],
        [ 'dms_env', 'DMS', 0.1, 2 ],
        [ 'dms_ha_h1', 'DMS', 0.1, 2 ],
        [ 'dms_ha_h3', 'DMS', 0.1, 2 ],
        [ 'dms_infa', 'DMS_min', 0.98, 1 ],
        [ 'dms_mapk1', 'DMS_SCH', 2.5, 1 ],
        [ 'dms_p53', 'DMS_null_etoposide', 1, 2 ],
        [ 'dms_pafa', 'DMS_kcat_km', 2300000, 1 ],
    ]

    for namespace, dms_name, cutoff, k_cutoff in settings:
        df = pd.read_csv(f'data/dms/{namespace}_esm.csv', delimiter=',')
        plt.figure()
        joint_grid = sns.jointplot(
            data=df,
            x='ESM_summed',
            y=dms_name,
            alpha=0.1,
            marginal_kws={ 'bins': 100, },
        )
        sns.scatterplot(
            x=df['ESM_summed'][df['ESM_vote'] >= k_cutoff],
            y=df[dms_name][df['ESM_vote'] >= k_cutoff],
            alpha=1.,
            color='#c41e3a',
            ax=joint_grid.ax_joint,
        )
        joint_grid.ax_joint.axhline(
            y=cutoff,
            color='gray',
            linestyle='--',
            zorder=10,
        )
        plt.savefig(f'figures/esm_scatter/{namespace}.png', dpi=500)
        plt.close()
