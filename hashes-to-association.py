#! /usr/bin/env python
"""
Given a list of sourmash signatures, find the common hashes and then
output a hash co-occurrence matrix that can be used to cluster samples.
"""
from __future__ import print_function

import argparse
import collections
import sourmash_lib.signature
import numpy


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
    args = p.parse_args()

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

    empty = set()
    for n, hashval1 in enumerate(hashlist):
        a = hashes_by_sig.get(hashval1, empty)
        for o, hashval2 in enumerate(hashlist):
            if n <= o:
                b = hashes_by_sig.get(hashval2, empty)
                common = len(a.intersection(b))
                pa[n][o] = common
                pa[o][n] = common

    pa /= float(n)

    print('\ndone! saving to:', args.output_name)

    with open(args.output_name, 'wb') as fp:
        numpy.save(fp, pa)

    with open(args.output_name + '.labels.txt', 'w') as fp:
        fp.write("\n".join(map(str, hashlist)))
            

if __name__ == '__main__':
    main()
