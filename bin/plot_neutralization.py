import matplotlib.pyplot as plt
import neutcurve
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import seaborn as sns
import sys

if __name__ == '__main__':
    fname = sys.argv[1]
    namespace = fname.split('/')[-1].split('.')[0]
    
    data, n_seen = [], {}
    with open(fname) as f:
        header = f.readline().rstrip().split('\t')
        conc_unit = header[1]
        for line in f:
            fields = line.rstrip().split('\t')
            virus = fields[0]
            concentration = float(fields[1])
            for condition, value in zip(header[2:], fields[2:]):
                #if condition not in { 's309', 'sotro', 'lc a1', 'a3a1', 'a6a1' }:
                #    continue
                if (virus, concentration, condition) not in n_seen:
                    n_seen[(virus, concentration, condition)] = 0
                n_seen[(virus, concentration, condition)] += 1
                data.append([
                    virus,
                    concentration,
                    condition,
                    n_seen[(virus, concentration, condition)],
                    float(value)/100.,
                ])

    df = pd.DataFrame(data, columns=[
        'virus',
        'concentration',
        'serum',
        'replicate',
        'fraction infectivity'
    ])

    fits = neutcurve.CurveFits(
        df,
        fixtop=False,
        fixbottom=0.,
    )

    print(fits.fitParams(average_only=False,))

    plt.figure()
    fig, axes = fits.plotViruses(
        xlabel=f'Concentration ({conc_unit})',
        legendfontsize=14,
        max_sera_per_subplot=7,
    )
    fig.set_size_inches(4.6, 4.7)
    plt.tight_layout()
    plt.savefig(f'figures/neutralization/{namespace}.svg')
    plt.close()

    exit()

    concentration_sorted = sorted(set(df['concentration']))
    conditions = sorted(set(df['condition']))

    plt.figure(figsize=(8, 5))
    sns.pointplot(
        data=df,
        x='concentration',
        y='value',
        hue='condition',
        col='virus',
        capsize=0.2,
    )
    #for condition in conditions:
    #    df_cond = df[df['condition'] == condition]
    #    popt, pcov = curve_fit(
    #        logifunc, df_cond['concentration'], df_cond['value'],
    #        p0=[ 1, 1, 10, 0]
    #    )
    #    plt.plot(
    #        df_cond['concentration'],
    #        logifunc(df_cond['concentration'], *popt),
    #        label=condition,
    #    )
    #    plt.scatter(
    #        df_cond['concentration'],
    #        df_cond['value'],
    #        label=condition,
    #    )
    #    #sns.regplot(
    #    #    data=df[df['condition'] == condition],
    #    #    x='concentration',
    #    #    y='value',
    #    #    logistic=True,
    #    #)
    #plt.xscale('log', base=10)
    plt.xlabel(f'Concentration ({conc_unit})')
    plt.ylabel('% Neutralization')
    plt.savefig(f'figures/neutralization/{namespace}.svg')
    plt.close()
