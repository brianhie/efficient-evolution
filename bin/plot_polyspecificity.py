import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def load_data():
    well_to_name = {}
    with open('data/polyspecificity/well_map.txt') as f:
        for line in f:
            well, name = line.rstrip().split('\t')
            well_to_name[well] = name

    data = []
    with open('data/polyspecificity/mfi_raw_table.txt') as f:
        for line in f:
            fields = line.rstrip().split('\t')

            if line.startswith('Plot '):
                well = fields[0].split()[-1]
                name = well_to_name[well]
                if well == 'A01':
                    columns = [ 'well', 'sample_name' ] + fields[1:]

            elif line.startswith('R2'):
                values = [ well, name ]
                for field in fields[1:]:
                    field = field.replace(',', '')
                    if field.endswith('%'):
                        field = float(field[:-1]) / 100.
                    values.append(float(field))
                data.append(values)

    df = pd.DataFrame(data, columns=columns)
    return df

if __name__ == '__main__':
    df = load_data()

    with open('data/polyspecificity/name_order.txt') as f:
        name_order = f.read().rstrip().split('\n')

    plt.figure(figsize=(15, 4))
    sns.barplot(data=df, x='sample_name', y='Median APC-A',
                order=name_order)
    sns.stripplot(
        x='sample_name',
        y='Median APC-A',
        data=df,
        order=name_order,
        color='black',
        size=6,
    )
    plt.yscale('log')
    plt.xticks(rotation=60)
    plt.tight_layout()
    plt.savefig('figures/polyspecificity.svg')
    plt.close()

    df.to_csv('data/polyspecificity/data.csv', sep=',')
