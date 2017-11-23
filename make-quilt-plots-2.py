#! /usr/bin/env python
"""
Given a list of sourmash signatures, find the common hashes and then
output a hash co-occurrence matrix that can be used to cluster samples.
"""
from __future__ import print_function

import random
import argparse
import collections
import sourmash_lib.signature
import numpy
import math
import sys
import os

import numpy
import scipy
import pylab
import scipy.cluster.hierarchy as sch


def mutinfo(total_n, n_common, n_a, n_b):
    mi = 0.0
    norm = True

    if n_common:
        mi = math.log(float(n_common) / n_a / n_b * total_n, 2)
        norm = -math.log(n_common / total_n, 2)
        if norm:
            mi /= norm
        else:
            mi = 1.0

    return mi #*n_common


def test_mutinfo():
    N = numpy.zeros((2, 2), dtype=numpy.float)
    N[0,0] = 1
    N[1,1] = 1
    N[0,1] = 1
    N[1,0] = 1

    common = 2
    n_a = 2
    n_b = 2
    total_n = 2

    assert mutinfo(total_n, common, n_a, n_b) == 1.0


def test_mutinfo_2():
    N = numpy.zeros((2, 2), dtype=numpy.float)
    N[0,0] = 1
    N[1,1] = 1
    N[0,1] = 0
    N[1,0] = 0

    common = 1
    n_a = 1
    n_b = 1
    total_n = 2

    assert mutinfo(total_n, common, n_a, n_b) == 0.0


def test_mutinfo_3():
    N = numpy.zeros((3, 3), dtype=numpy.float)
    N[0,0] = 1
    N[1,1] = 1
    N[0,1] = 0
    N[1,0] = 0
    N[2,2] = 1

    pmi(N, 0, 0)

    common = 1
    n_a = 1
    n_b = 1
    total_n = 2

    assert mutinfo(total_n, common, n_a, n_b) == 0.0


def plot_matrix(D):
    """Build a composite plot showing dendrogram + distance matrix/heatmap.

    Returns a matplotlib figure."""
    fig = pylab.figure(figsize=(7, 7))

    # compute dendrogram (but don't plot it)
    Y = sch.linkage(D, method='single')  # centroid

    Z1 = sch.dendrogram(Y, orientation='left', no_labels=True, no_plot=True)

    # plot matrix
    axmatrix = fig.add_axes([0.1, 0.1, 0.75, 0.75])

    # (this reorders D by the clustering in Z1)
    idx1 = Z1['leaves']
    D = D[idx1, :]
    D = D[:, idx1]

    # show matrix
    im = axmatrix.matshow(D, aspect='auto', origin='lower',
                          cmap=pylab.cm.YlGnBu, vmax=1.0, vmin=0.0)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([0.852, 0.1, 0.02, 0.75])
    pylab.colorbar(im, cax=axcolor)

    return fig


def main():
    p = argparse.ArgumentParser()
    p.add_argument('inp_signatures', nargs='+')
    p.add_argument('-k', '--ksize', type=int, default=31)
    p.add_argument('--scaled', type=int, default=100000)
    p.add_argument('--query', nargs='+',
                   help='produce quilt plots for each query sig')
    args = p.parse_args()

    total_n = len(args.inp_signatures)
    counts = collections.Counter()

    if not args.query:
        print('must specify query')
        sys.exit(-1)

    print('loading signatures from', len(args.inp_signatures), 'files')
    sig_hashes = {}
    hashes_by_sig = collections.defaultdict(set)
    all_hashes = set()
    for n, filename in enumerate(args.inp_signatures):
        print('... {}'.format(n + 1), end='\r')
        sig = sourmash_lib.load_one_signature(filename, ksize=args.ksize)
        mh = sig.minhash.downsample_scaled(args.scaled)
        hashes = mh.get_mins()
        all_hashes.update(hashes)

        sig_hashes[filename] = hashes

        for hashval in hashes:
            hashes_by_sig[hashval].add(filename)

    print('\n...done. Now calculating associations.')

    ## now, for each query...
    listfp = open('list.txt', 'wt')
    listassocfp = open('list.assoc.txt', 'wt')
    for filename in args.query:
        #print('...loading query {}'.format(filename))
        sig = sourmash_lib.signature.load_one_signature(filename,
                                                        ksize=args.ksize)
        mh = sig.minhash.downsample_scaled(args.scaled)
        query_hashes = set(mh.get_mins())

        # intersect with the database of samples
        query_hashes.intersection_update(all_hashes)

        if len(query_hashes) <= 1:
            print('SKIPPING {}, no intersect hashes'.format(filename))
            continue

        # build matrix
        pa = numpy.zeros((len(query_hashes), len(query_hashes)),
                          dtype=numpy.float)

        print('calculating matrix {} x {}'.format(len(query_hashes),
                                                  len(query_hashes)))

        hashlist = list(query_hashes)

        for n, hashval1 in enumerate(hashlist):
            if n % 1000 == 0 and n:
                print('...', n)
            a = hashes_by_sig.get(hashval1)
            for o, hashval2 in enumerate(hashlist):
                if o <= n:
                    b = hashes_by_sig.get(hashval2)
                    common = len(a.intersection(b))

                    mi = mutinfo(total_n, common, len(a), len(b))
                    pa[n][o] = mi
                    pa[o][n] = mi

        print('matrix scale: min {}, max {}'.format(pa.min(), pa.max()))

        output_name = os.path.basename(filename) + '.quilt'

        with open(output_name, 'wb') as fp:
            numpy.save(fp, pa)

        with open(output_name + '.labels.txt', 'w') as fp:
            fp.write("\n".join(map(str, hashlist)))

        fig = plot_matrix(pa)
        fig.savefig(output_name + '.pdf')
        pylab.close()

        ## calculate a single number representing ...something.
        pamin, pamax = pa.min(), pa.max()
        pa -= pa.min()
        if pa.max() == 0.0:
            pa = numpy.ones(pa.shape, dtype=numpy.float)
        else:
            pa /= pa.max()

        pa_zero_val = numpy.sum(numpy.square(pa)) / len(query_hashes)**2
        print(os.path.basename(filename), pa_zero_val, pamin, pamax, file=listfp)

        ## @CTB look at associations
        assoc_hashes = set()
        for h in all_hashes:
            if h in query_hashes:
                assoc_hashes.add(h)
                continue
            a = hashes_by_sig.get(h)
            for qh in query_hashes:
                b = hashes_by_sig.get(qh)
                common = len(a.intersection(b))

                mi = mutinfo(total_n, common, len(a), len(b))
                if mi > 0.5:
                    assoc_hashes.add(h)

        print('got {} assoc hashes'.format(len(assoc_hashes)))
        print('downsampling to 500')

        hashlist = random.sample(assoc_hashes, 500)
        hashlist = list(hashlist)

        pa = numpy.zeros((len(hashlist), len(hashlist)),
                          dtype=numpy.float)

        for n, hashval1 in enumerate(hashlist):
            if n % 1000 == 0 and n:
                print('...', n)
            a = hashes_by_sig.get(hashval1)
            for o, hashval2 in enumerate(hashlist):
                if o <= n:
                    b = hashes_by_sig.get(hashval2)
                    common = len(a.intersection(b))

                    mi = mutinfo(total_n, common, len(a), len(b))
                    pa[n][o] = mi
                    pa[o][n] = mi

        print('matrix scale: min {}, max {}'.format(pa.min(), pa.max()))

        output_name = os.path.basename(filename) + '.assoc.quilt'

        with open(output_name, 'wb') as fp:
            numpy.save(fp, pa)

        with open(output_name + '.labels.txt', 'w') as fp:
            fp.write("\n".join(map(str, hashlist)))

        fig = plot_matrix(pa)
        fig.savefig(output_name + '.pdf')
        pylab.close()
        
        ## calculate a single number representing ...something.
        pamin, pamax = pa.min(), pa.max()
        pa -= pa.min()
        if pa.max() == 0.0:
            pa = numpy.ones(pa.shape, dtype=numpy.float)
        else:
            pa /= pa.max()

        pa_zero_val = numpy.sum(numpy.square(pa)) / len(hashlist)**2
        print(os.path.basename(filename), pa_zero_val, pamin, pamax, file=listassocfp)
            

if __name__ == '__main__':
    main()
