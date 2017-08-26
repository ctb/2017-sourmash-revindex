#! /usr/bin/env python
import argparse
import collections
import sys
sys.path.insert(0, 'kraken/')
sys.path.insert(0, 'revindex/')

import revindex_utils
import lca_json

SCALED=100000


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('revindex')
    p.add_argument('lca_json')
    args = p.parse_args()

    idx = revindex_utils.HashvalRevindex(args.revindex)
    lca_db = lca_json.LCA_Database(args.lca_json)
    taxfoo, hashval_to_lca, scaled = lca_db.get_database(args.ksize, SCALED)
    
    cnt = collections.Counter()
    for k, v in idx.hashval_to_sigids.items():
        cnt[k] += len(v)

    total = 0
    found = 0
    unknown = collections.defaultdict(int)
    for k, c in cnt.most_common():
        if c < 10:              # break when we hit things in < 10 samples.
            break
        total += 1
        lca = hashval_to_lca.get(k)
        if lca:
            print('hash {}, in {} samples; lineage: {}'.format(k, c, ";".join(taxfoo.get_lineage(lca))))
            found += 1
        else:
            unknown[c] += 1

    print('found', found, 'of', total)
    print('distribution of unknowns:')
    print('commonality n')

    sofar = 0
    for k, cnt in sorted(unknown.items()):
        sofar += cnt
        print('{} {} {}'.format(k, cnt, sofar))


if __name__ == '__main__':
    main()
