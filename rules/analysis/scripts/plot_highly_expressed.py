#!/usr/bin env python

from collections import OrderedDict
import sys
import argparse
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import yaml
import pandas as pd
import numpy as np



def higly_expressed_yaml(E, F, ntop=10):

    order = E.values.sum(1).argsort()[::-1]
    keep = order[:ntop]
    E.f = E.iloc[keep,]
    F.f = F.loc[E.f.index,:]
    df = pd.concat([E.f, F.f[['gene_name', 'gene_biotype']]], axis=1)
    
    

    cats = []
    keys = OrderedDict()
    default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c',
                      '#8085e9','#f15c80', '#e4d354', '#2b908f',
                      '#f45b5b', '#91e8e1']
    for i, k in enumerate(df.index):
        name = df.loc[k, 'gene_name'] + ': ' + df.loc[k, 'gene_biotype']
        color = default_colors[i]
        keys[k] = {'color': color, 'name': name}
    
    # Config for the plot
    pconfig = {
        'id': 'gene_high',
        'title': 'Higly Expressed Genes',
        'ylab': 'TPM'
    }
    
    section = {}
    section['id'] = 'gene_high'
    section['section_name'] = 'Higly Expressed Genes'
    section['description'] = 'Summary of the 10 most highly expressed genes.'
    section['plot_type'] = 'bargraph'
    
    section['pconfig'] = pconfig
    section['categories'] = keys
    data_dict = OrderedDict()
    for k, v in df.to_dict().items():
        if not k in ['gene_name', 'gene_biotype']:
            data_dict[k] = v
    section['data'] = [data_dict]

    return section


def argparser():
    parser = argparse.ArgumentParser(description="Higly expressed genes figure for QC report")
    parser.add_argument("exprs")
    parser.add_argument(
        "--sample-info",
        help="Optional sample sheet. Will subset expr table if needed",
        dest="samples",
    )
    parser.add_argument(
        "--feature-info",
        help="Required feature info. Will subset expr table if needed",
        dest="features", required=True,
    )
    parser.add_argument(
        "-o ",
        "--output",
        default="biotypes_mqc.pyaml",
        help="Output filename. Will default to biotypes_mqc.yaml",
    )

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = argparser()

    E = pd.read_csv(args.exprs, sep="\t", index_col=0)
    E.columns = E.columns.astype(str)
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        S.index = S.index.astype(str)
        if not E.columns.isin(S.index).all():
            raise ValueError("missing samples in sample info!")
        S = S.loc[E.columns, :]
    
    F = pd.read_csv(args.features, sep="\t", index_col=0)
    if not E.index.isin(F.index).all():
        warnings.warn("missing annotations in feature info!")
    F = F.loc[E.index, :]

    if not 'gene_biotype' in F.columns or 'gene_name' not in F.columns:
        print(F.columns)
        raise ValueError('Feature info needs columns `gene_biotype` and `gene_name` !')


    section =  higly_expressed_yaml(E, F)

    with open(args.output, 'w') as fh:
        yaml.dump(section, fh, default_flow_style=False)
