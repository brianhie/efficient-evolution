import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as ss
import seaborn as sns

if __name__ == '__main__':
    df = pd.read_csv('data/ic50_kd_relationship.txt', sep='\t')

    df = df[df['Design'] != 'WT']
    antibodies = sorted(set(df['Antibody']))

    print(ss.spearmanr(df['Kd Fold change'], df['IC50 Fold Change']))

    fig = plt.figure()
    sns.scatterplot(
        data=df,
        x='Kd Fold change',
        y='IC50 Fold Change',
        hue='Antibody',
        s=200,
    )
    plt.xscale('log')
    plt.yscale('log')
    fig.set_size_inches(4, 4)
    plt.savefig('figures/ic50_kd.svg')
    plt.close()
    
    fig = plt.figure()
    sns.scatterplot(
        data=df[(df['Antibody'] != 'mAb114 UCA') &
                (df['Antibody'] != 'C143')],
        x='Kd Fold change',
        y='IC50 Fold Change',
        hue='Antibody',
        s=200,
    )
    plt.xscale('log')
    plt.yscale('log')
    fig.set_size_inches(4, 4)
    plt.savefig('figures/ic50_kd_subset.svg')
    plt.close()
    
    fig, axes = plt.subplots(1, 3)
    for ab_idx, antibody in enumerate(antibodies):
        if antibody == 'mAb114 UCA':
            continue
        sns.scatterplot(
            data=df[df['Antibody'] == antibody],
            x='Kd Fold change',
            y='IC50 Fold Change',
            ax=axes[ab_idx]
        )
    fig.set_size_inches(10, 5)
    plt.savefig('figures/ic50_kd_subplots.svg')
    plt.close()
