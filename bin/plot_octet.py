import glob
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys

def load_file(fname, namespace=None):
    if namespace is None:
        namespace = fname.split('/')[-1].split('.')[0]

    data = []
    with open(fname) as f:
        for line in f:
            if line.startswith('Vers') or \
               line.startswith('Conc') or \
               line.startswith('Start') or \
               line.startswith('Stop') or \
               line.startswith('Time'):
                continue
        
            fields = line.rstrip().split(',')
            time = float(fields[0])
            value = float(fields[1])
            data.append([ namespace, time, value ])

    return data


if __name__ == '__main__':
    dirname = sys.argv[1]

    data = []
    for fname in glob.glob(f'{dirname}/*.csv'):
        if fname.endswith('kineticanalysistableresults.csv'):
            continue
        if '/c055-' in fname:
            continue
        print(fname)
        data += load_file(fname)

    df = pd.DataFrame(data, columns=[ 'namespace', 'time', 'value' ])

    #new_data = []
    #for namespace in set(df['namespace']):
    #    df_subset = df[df['namespace'] == namespace]
    #    max_val = max(df_subset['value'])
    #    for namespace, time, value in zip(df_subset.namespace,
    #                                      df_subset.time,
    #                                      df_subset.value):
    #        new_data.append([ namespace, time, value / max_val ])
    #df = pd.DataFrame(new_data, columns=[ 'namespace', 'time', 'value' ])

    #df = df[df.time > 2500]

    plt.figure(figsize=(5, 5))
    sns.lineplot(data=df, x='time', y='value', hue='namespace',
                 linewidth=3.)
    plt.legend(loc='lower right')
    plt.savefig('figures/plot_octet.svg')
    plt.close()
