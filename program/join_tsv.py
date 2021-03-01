import click
from sys import stdout
import pandas as pd


@click.command()
@click.argument('tsv', nargs=-1)
@click.option('-o', '--out', help='Output joined tsv')
def join_tsv(tsv, out):
    df_list = []
    for file in tsv:
        df = pd.read_csv(file, sep='\t', index_col=0)
        df_list.append(df)

    joined_df = pd.concat(df_list, axis=1, sort=False)

    out_file = out if out else stdout

    joined_df.index.name = 'Assemblies'
    joined_df.reset_index(level=0, inplace=True)

    # joined_df['Assemblies'] = joined_df.index

    joined_df.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    join_tsv()
