#! /usr/bin/env python
import argparse
from revindex_utils import HashvalRevindex
import collections

def main():
    p = argparse.ArgumentParser()
    p.add_argument('loadname')
    args = p.parse_args()

    idx = HashvalRevindex(args.loadname)

    cnt = collections.Counter()
    for k, v in idx.hashval_to_sigids.items():
        cnt[k] += len(v)

    for k, c in cnt.most_common(100):
        print(k, c)


if __name__ == '__main__':
    main()
