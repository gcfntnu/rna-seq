#!/usr/bin env python

import sys
import warnings

warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import argparse
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt


def var_filter(E, f=0.5):
    E = E.loc[E.sum(1) > 10, :]
    o = E.var(axis=1).argsort()[::-1]
    n = int(np.floor(E.shape[0] * f))
    E = E.iloc[o[:n], :]
    return E


def pca(E, max_comp=20):
    x = E.values
    max_comp = min(max_comp, min(x.shape))
    x = x - x.mean(0)
    u, s, vt = np.linalg.svd(x, 0)
    T = u * s
    P = vt.T
    T = pd.DataFrame(
        T[:, :max_comp],
        index=E.index,
        columns=["PC_{}".format(i + 1) for i in range(max_comp)],
    )
    P = pd.DataFrame(
        P[:, :max_comp],
        index=E.columns,
        columns=["PC_{}".format(i + 1) for i in range(max_comp)],
    )
    return T, P


def pairs_scores(T, S=None, n_comp=4):
    n_comp = min(T.shape[1], n_comp)
    T = T.iloc[:, :n_comp]
    df = pd.concat([T, S], axis="columns")
    sns.set(style="ticks")
    hue = "Sample_Group" if "Sample_Group" in df.columns else None
    ax = sns.pairplot(df, vars=T.columns, hue=hue)
    return ax


def argparser():
    parser = argparse.ArgumentParser(description="PCA figure for QC report")
    parser.add_argument("exprs")
    parser.add_argument(
        "--sample-info",
        help="Optional sample sheet. Will subset expr table if needed",
        dest="samples",
    )
    parser.add_argument(
        "--feature-info",
        help="Optional sample sheet. Will subset expr table if needed",
        dest="features",
    )
    parser.add_argument(
        "-o ",
        "--output",
        default="pca_mqc.png",
        help="Output filename. Will default to pca.png",
    )
    parser.add_argument("--include-comp-test", action="store_true", dest="test_comp")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argparser()

    E = pd.read_csv(args.exprs, sep="\t", index_col=0)
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        if not E.columns.isin(S.index).all():
            raise ValueError("missing samples in sample info!")
        S = S.loc[E.columns, :]
    if args.features is not None:
        F = pd.read_csv(args.features, sep="\t", index_col=0)
        if not E.index.isin(F.index).all():
            warnings.warn("missing annotations in feature info!")
        F = F.loc[E.index, :]

    F = var_filter(E)
    T, P = pca(F.T)
    p = pairs_scores(T, S)
    plt.subplots_adjust(top=0.9)
    p.fig.suptitle("PCA on variance transformed gene expression data")
    p.savefig(args.output)
