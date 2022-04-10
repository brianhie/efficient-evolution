import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

if __name__ == '__main__':
    df = pd.read_csv('data/dms/enrichments.txt', sep='\t')

    data = []
    for protein, hit_rate, null_pct in zip(
            df['Protein'],
            df['Hit rate'],
            df['Null percentage'],
    ):
        data.append([ protein, 'Hit rate', hit_rate ])
        data.append([ protein, 'Null percentage', null_pct ])
        
    df = pd.DataFrame(data, columns=[ 'protein', 'pct_type', 'pct' ])

    plt.figure()
    ax = sns.catplot(
        data=df,
        kind='bar',
        x='protein',
        y='pct',
        hue='pct_type',
        aspect=1.5,
    )
    plt.xlabel('Protein')
    plt.ylabel('%')
    plt.savefig('figures/dms_enrichments.svg')
    plt.close()
