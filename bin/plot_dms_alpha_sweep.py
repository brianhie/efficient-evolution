import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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
        [ 'dms_adrb2', 'DMS_0.625', 2.8 ],
        [ 'dms_bla', 'DMS_amp_2500_(b)', 0.01 ],
        [ 'dms_env', 'DMS', 0.1 ],
        [ 'dms_ha_h1', 'DMS', 0.1 ],
        [ 'dms_ha_h3', 'DMS', 0.1 ],
        [ 'dms_infa', 'DMS_min', 0.98 ],
        [ 'dms_mapk1', 'DMS_SCH', 2.5 ],
        [ 'dms_p53', 'DMS_null_etoposide', 1 ],
        [ 'dms_pafa', 'DMS_kcat_km', 2300000 ],
    ]

    data = []

    for alpha in [ 0, 0.1, 0.2, 0.5, 0.7, 0.9,]:
    
        dms_lm_fname = f'target/log/dms_alpha{alpha}.log'
        df_lm = load_dms_lm(dms_lm_fname)
    
        for namespace, dms_name, cutoff in settings:
            fname = f'data/dms/{namespace}.csv'
            with open(fname):
                df_dms = pd.read_csv(fname, delimiter=',')
            fitness_finite = df_dms[dms_name][np.isfinite(df_dms[dms_name])]
            fitness_max = max(fitness_finite)
            fitness_99pct = np.percentile(fitness_finite, 99)
            fitness_95pct = np.percentile(fitness_finite, 95)
            fitness_90pct = np.percentile(fitness_finite, 90)
            fitness_75pct = np.percentile(fitness_finite, 75)
            fitness_min = min(fitness_finite)

            df_setting = df_lm[
                (df_lm['namespace'] == namespace) &
                (df_lm['dms_name'] == dms_name)
            ]
            
            n_positive = len(df_setting[
                (df_setting['fitness'] > cutoff)
            ])
            n_total = len(df_setting)
            if n_total == 0:
                continue
            lm_guided_pct = n_positive / n_total * 100.

            fitness_raw = np.array(df_setting['fitness'])
            fitness_norm = pct_norm(fitness_raw, fitness_min, fitness_max)
            fitness_norm99 = pct_norm(fitness_raw, fitness_min, fitness_99pct)
            fitness_norm95 = pct_norm(fitness_raw, fitness_min, fitness_95pct)
            fitness_norm90 = pct_norm(fitness_raw, fitness_min, fitness_90pct)
            fitness_norm75 = pct_norm(fitness_raw, fitness_min, fitness_75pct)
            fitness_norm99[fitness_norm99 > 1] = 1
            fitness_norm95[fitness_norm95 > 1] = 1
            fitness_norm90[fitness_norm90 > 1] = 1
            fitness_norm75[fitness_norm75 > 1] = 1
                
            data.append([
                namespace,
                alpha,
                n_positive,
                n_total,
                lm_guided_pct,
                max(fitness_norm),
                max(fitness_norm99),
                max(fitness_norm95),
                max(fitness_norm90),
                max(fitness_norm75),
            ])

    dms_lm_fname = 'target/log/dms.log'
    df_lm = load_dms_lm(dms_lm_fname)
    max_models = 3

    for namespace, dms_name, cutoff in settings:
        df_setting = df_lm[
            (df_lm['namespace'] == namespace) &
            (df_lm['dms_name'] == dms_name)
        ]
        
        fname = f'data/dms/{namespace}.csv'
        with open(fname):
            df_dms = pd.read_csv(fname, delimiter=',')
        fitness_finite = df_dms[dms_name][np.isfinite(df_dms[dms_name])]
        fitness_max = max(fitness_finite)
        fitness_99pct = np.percentile(fitness_finite, 99)
        fitness_95pct = np.percentile(fitness_finite, 95)
        fitness_90pct = np.percentile(fitness_finite, 90)
        fitness_75pct = np.percentile(fitness_finite, 75)
        fitness_min = min(fitness_finite)

        for k in range(max_models):
            n_positive = len(df_setting[
                (df_setting['n_models'] > k) &
                (df_setting['fitness'] > cutoff)
            ])
            n_total = len(df_setting[
                (df_setting['n_models'] > k)
            ])
            if n_total == 0:
                continue
            lm_guided_pct = n_positive / n_total * 100.

            fitness_raw = np.array(df_setting['fitness'])
            fitness_norm = pct_norm(fitness_raw, fitness_min, fitness_max)
            fitness_norm99 = pct_norm(fitness_raw, fitness_min, fitness_99pct)
            fitness_norm95 = pct_norm(fitness_raw, fitness_min, fitness_95pct)
            fitness_norm90 = pct_norm(fitness_raw, fitness_min, fitness_90pct)
            fitness_norm75 = pct_norm(fitness_raw, fitness_min, fitness_75pct)
            fitness_norm99[fitness_norm99 > 1] = 1
            fitness_norm95[fitness_norm95 > 1] = 1
            fitness_norm90[fitness_norm90 > 1] = 1
            fitness_norm75[fitness_norm75 > 1] = 1
                
            data.append([
                namespace,
                f'n_models{k+1}',
                n_positive,
                n_total,
                lm_guided_pct,
                max(fitness_norm),
                max(fitness_norm99),
                max(fitness_norm95),
                max(fitness_norm90),
                max(fitness_norm75),
            ])
            
    df = pd.DataFrame(data, columns=[
        'namespace',
        'alpha',
        'n_positive',
        'n_total',
        'value',
        'fitness_max',
        'fitness_99pct',
        'fitness_95pct',
        'fitness_90pct',
        'fitness_75pct',
    ])
    n_bars = len(set(df['alpha']))

    val_names = [
        'n_total',
        'value',
        'fitness_max',
        'fitness_99pct',
    ]
    palettes = [
        sns.cubehelix_palette(n_colors=n_bars, start=2.8, rot=0.1),
        sns.cubehelix_palette(n_colors=n_bars),
        sns.cubehelix_palette(n_colors=n_bars, rot=-0.4),
        sns.cubehelix_palette(n_colors=n_bars, rot=-0.6),
    ]

    for val_name, palette in zip(val_names, palettes):
        plt.figure(figsize=(10, 3 + 2*(val_name == 'value')))
        sns.barplot(data=df, x='namespace', y=val_name, hue='alpha',
                    palette=palette)
        if val_name in { 'n_total' }:
            plt.yscale('log')
        plt.savefig(f'figures/dms_alpha_sweep_{val_name}.svg')
        plt.close()
