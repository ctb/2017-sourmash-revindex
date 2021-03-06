#! /usr/bin/env python
"""
Given a list of sourmash signatures, find the common hashes and then
output a hash co-occurrence matrix that can be used to cluster samples.
"""
from __future__ import print_function

import sys
import argparse
import collections
import numpy

import sourmash_lib.signature

def main():
    p = argparse.ArgumentParser()
    p.add_argument('inp_signatures', nargs='+')
    p.add_argument('-k', '--ksize', type=int, default=31)
    p.add_argument('--scaled', type=int, default=100000)
    p.add_argument('-o', '--output-name')
    p.add_argument('--threshold', type=int, default=2)
    p.add_argument('--max-threshold', type=int, default=None)
    p.add_argument('--intersect', nargs='+',
                   help='only use hashes in the given files')
    p.add_argument('--query', nargs='+')
    p.add_argument('--query-association', type=float, default=0.8)
    args = p.parse_args()

    assert args.query

    counts = collections.Counter()

    intersect_hashes = set()
    if args.intersect:
        for n, filename in enumerate(args.intersect):
            print('...loading intersect {}'.format(n + 1), end='\r')
            sig = sourmash_lib.signature.load_one_signature(filename,
                                                            ksize=args.ksize)
            mh = sig.minhash.downsample_scaled(args.scaled)
            hashes = mh.get_mins()
            intersect_hashes.update(hashes)
        print('')

    print('loading signatures from', len(args.inp_signatures), 'files')
    sig_hashes = {}
    hashes_by_sig = collections.defaultdict(set)
    for n, filename in enumerate(args.inp_signatures):
        print('... {}'.format(n + 1), end='\r')
        sig = sourmash_lib.load_one_signature(filename, ksize=args.ksize)
        mh = sig.minhash.downsample_scaled(args.scaled)
        hashes = mh.get_mins()

        if intersect_hashes:
            hashes = set(hashes)
            hashes.intersection_update(intersect_hashes)

        sig_hashes[filename] = hashes

        for hashval in hashes:
            hashes_by_sig[hashval].add(filename)
            counts[hashval] += 1

    print('\n...done. Now calculating associations.'.format(args.threshold))

    n = 0
    abundant_hashes = set()
    for hashval, count in counts.most_common():
        if args.max_threshold and count > args.max_threshold:
            continue

        if count < args.threshold:
            break

        n += 1
        abundant_hashes.add(hashval)

    # filter by association in a set of queries - first, load queries
    print('filtering associations by query from {}'.format(args.query))
    query_hashes = set()
    for query_filename in args.query:
        query_sig = sourmash_lib.load_one_signature(query_filename,
                                                    ksize=args.ksize)
        query_mh = query_sig.minhash.downsample_scaled(args.scaled)
        query_hashes.update(query_mh.get_mins())

    query_hashes.intersection_update(abundant_hashes)

    # second, filter out only those with some assoc!
    new_hashes = set()
    for qhash in query_hashes:
        a = hashes_by_sig[qhash]
        for hashval in abundant_hashes:
            b = hashes_by_sig[hashval]

            common = len(a.intersection(b))
            frac = 2*common / float(len(a) + len(b))

            if frac > args.query_association:
                new_hashes.add(hashval)

    print('filtered hashes from {} to {}'.format(len(abundant_hashes), len(new_hashes)))
    before = len(new_hashes)
    new_hashes -= query_hashes
    after = len(new_hashes)
    print('now down to {} (-{})'.format(len(new_hashes), before - after))

    abundant_hashes = new_hashes

    if not args.output_name:
        sys.exit(0)

    # sort for no particular reason
    hashlist = list(sorted(abundant_hashes))
    hashdict = {}
    for n, hashval in enumerate(hashlist):
        hashdict[hashval] = n                   # hash -> index in hashlist

    print('found', n, 'hashes from', len(args.inp_signatures), 'signatures')
    print('min threshold: {}'.format(args.threshold))

    pa = numpy.zeros((len(abundant_hashes), len(abundant_hashes)),
                      dtype=numpy.float)

    print('calculating matrix {} x {}'.format(len(abundant_hashes),
                                              len(abundant_hashes)))

    for n, hashval1 in enumerate(hashlist):
        if n % 1000 == 0:
            print('...', n)
        a = hashes_by_sig.get(hashval1)
        for o, hashval2 in enumerate(hashlist):
            if o <= n:
                b = hashes_by_sig.get(hashval2)
                common = len(a.intersection(b))
                frac = 2*common / float(len(a) + len(b))
                pa[n][o] = frac
                pa[o][n] = frac

    print('XXX', pa.min(), pa.max())
    print('\ndone! saving to:', args.output_name)

    if args.output_name:
        with open(args.output_name, 'wb') as fp:
            numpy.save(fp, pa)

        with open(args.output_name + '.labels.txt', 'w') as fp:
            fp.write("\n".join(map(str, hashlist)))
            

if __name__ == '__main__':
    main()
