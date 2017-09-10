#! /usr/bin/env python
"""
Given a list of sourmash signatures, find the common hashes and
then output a file x hash matrix that can be used to cluster samples
by common hashes.
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
    args = p.parse_args()

    counts = collections.Counter()

    print('loading signatures from', len(args.inp_signatures), 'files')
    for filename in args.inp_signatures:
        sig = sourmash_lib.signature.load_one_signature(filename,
                                                        select_ksize=args.ksize)
        mh = sig.minhash.downsample_scaled(args.scaled)
        hashes = mh.get_mins()

        for k in hashes:
            counts[k] += 1

    n = 0
    abundant_hashes = set()
    for hash, count in counts.most_common():
        if count >= args.threshold:
            n += 1
            abundant_hashes.add(hash)

    print('found', n, 'hashes from', len(args.inp_signatures), 'signatures')
    print('min threshold: {}'.format(args.threshold))

    # go over the files again, this time creating an n x n_files matrix
    # with 0 etc.
    pa = numpy.zeros((len(args.inp_signatures), len(abundant_hashes)),
                      dtype=numpy.int)

    # sort for no particular reason
    hashlist = list(sorted(abundant_hashes))
    hashdict = {}
    for n, k in enumerate(hashlist):
        hashdict[k] = n                   # hash -> index in hashlist
                         
    print('load x 2 signatures from', len(args.inp_signatures), 'files')
    for fn, filename in enumerate(args.inp_signatures):
        sig = sourmash_lib.signature.load_one_signature(filename,
                                                        select_ksize=args.ksize)
        mh = sig.minhash
        mh = mh.downsample_scaled(args.scaled)
        hashes = mh.get_mins()

        x = abundant_hashes.intersection(hashes)
        for hashval in x:
            idx = hashdict[hashval]
            pa[fn][idx] = 1

    print('saving to:', args.output_name)

    with open(args.output_name, 'wb') as fp:
        numpy.save(fp, pa)

    with open(args.output_name + '.labels.txt', 'w') as fp:
        fp.write("\n".join(args.inp_signatures))
            

if __name__ == '__main__':
    main()
