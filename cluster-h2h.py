#! /usr/bin/env python
import argparse
import numpy as np
import sklearn.cluster as cluster
import hdbscan
import numpy as np
import pandas as pd

import sourmash_lib.fig


SELECT_TOP_N=100


# create a data frame that connects sample_index, cluster label, hash value, and cluster size
def make_sample_df(cluster_labels, sample_labels):
    df_cl = pd.DataFrame(cluster_labels)
    df_cl.index.name = 'sample_index'
    df_cl.columns = ['cluster']
    df_cl['hashval'] = pd.DataFrame(list(map(int, sample_labels)))
    df_cl = df_cl[df_cl.cluster != -1]

    df_cl['cluster_size'] = df_cl.groupby('cluster')['cluster'].transform('count')

    return df_cl


def main():
    p = argparse.ArgumentParser()
    p.add_argument('h2n_matrix')
    args = p.parse_args()

    print('loading matrix {}'.format(args.h2n_matrix))
    mat, sample_labels = sourmash_lib.fig.load_matrix_and_labels(args.h2n_matrix)
    cluster_labels = hdbscan.HDBSCAN(min_cluster_size=20).fit_predict(mat)

    df = make_sample_df(cluster_labels, sample_labels)

    cluster_sizes = list(set(df.cluster_size))
    cluster_sizes.sort()
    cluster_sizes.reverse()
    cluster_sizes = cluster_sizes[:100]
    
    cluster_ids = list(df[df.cluster_size.isin(cluster_sizes)].cluster)

    for cluster_id in cluster_ids:
        mh = sourmash_lib.MinHash(0, 51, scaled=100000)
        mh.add_many(list(df[df.cluster == cluster_id].hashval))

        sig = sourmash_lib.SourmashSignature('', mh, name='cluster{}'.format(cluster_id))
        filename = '{}.cluster{}.sig'.format(args.h2n_matrix, cluster_id)
        sourmash_lib.save_signatures([sig], open(filename, 'wt'))


if __name__ == '__main__':
    main()
