#! /usr/bin/env python
"""
Given one or more LCA databases, and a revindex, apply LCA
classification to those hashes with abundance higher than
--abundance-threshold in more than --sample-threshold samples.
Then output statistics about how many remain unclassified.
"""
import argparse
import collections
import sys

import revindex_utils
from sourmash_lib.lca import lca_utils

SCALED=100000
DEFAULT_KSIZE=31
DEFAULT_SAMPLE_THRESHOLD=2
DEFAULT_ABUND_THRESHOLD=1


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    p.add_argument('--sample-threshold',
                   default=DEFAULT_SAMPLE_THRESHOLD, type=int)
    p.add_argument('--abundance-threshold',
                   default=DEFAULT_ABUND_THRESHOLD, type=int)
    p.add_argument('revindex')
    p.add_argument('db', nargs='+')
    args = p.parse_args()

    idx = revindex_utils.HashvalRevindex(args.revindex)

    lca_db_list, ksize, scaled = lca_utils.load_databases(args.db, SCALED)
    
    cnt = collections.Counter()
    for k, v in idx.hashval_to_abunds.items():
        cnt[k] += len([abund for abund in v \
                       if abund >= args.abundance_threshold])

    total = 0
    found = 0
    unknown = collections.defaultdict(int)
    for hashval, count in cnt.most_common():
        # break when we hit things in < 10 samples.        
        if count < args.sample_threshold:
            break
        total += 1
        lca_set = set()

        for lca_db in lca_db_list:
            lineages = lca_db.get_lineage_assignments(hashval)
            lca_set.update(lineages)

        if not lca_set:
            unknown[count] += 1
            continue

        assert lca_set, lca_set

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        tree = lca_utils.build_tree(lca_set)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = lca_utils.find_lca(tree)

        print('hash {}, in {} samples; lineage: {}'.format(hashval, count, ";".join(lca_utils.zip_lineage(lca))), file=sys.stderr)
        found += 1

    print('found', found, 'of', total, file=sys.stderr)
    print('outputting distribution of unknowns', file=sys.stderr)
    print('commonality,n,sum_n')

    sofar = 0
    for k, cnt in sorted(unknown.items()):
        sofar += cnt
        print('{},{},{}'.format(k, cnt, sofar))


if __name__ == '__main__':
    main()
