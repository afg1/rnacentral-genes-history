#!/usr/bin/env python

import polars as pl
import click


def transform_organism_name(name: str) -> str:
    """Transform organism name to match directory format"""
    import re
    transformed = name.lower()
    transformed = re.sub(r'[^a-z0-9]+', '_', transformed)
    transformed = transformed.strip('_')
    return transformed


@click.command()
@click.argument('input_csv', type=click.Path(exists=True))
@click.argument('output_csv', type=click.Path())
def main(input_csv, output_csv):
    df = pl.read_csv(input_csv)
    df = df.with_columns(
        pl.col('organism_name').map_elements(transform_organism_name).alias('transformed_name')
    )
    df.write_csv(output_csv)


if __name__ == '__main__':
    main()