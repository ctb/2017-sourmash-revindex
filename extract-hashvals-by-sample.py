#! /usr/bin/env python
"""
"""
import sys, os
import sourmash_lib, sourmash_lib.signature
import argparse
from pickle import dump, load
from collections import defaultdict

def traverse_find_sigs(dirnames):
    for dirname in dirnames:
        for root, dirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.sig') or name.endswith('.sbt'):
                    fullname = os.path.join(root, name)
                    yield fullname


def main():
    p = argparse.ArgumentParser()
    p.add_argument('savename')
    p.add_argument('sigs', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--scaled', default=10000, type=int)

    p.add_argument('--traverse-directory', action='store_true')

    args = p.parse_args()

    # track hashval -> list of signature IDs
    hashval_to_sigids = defaultdict(list)

    # track hashval -> list of abundances
    hashval_to_abunds = defaultdict(list)

    # track sigIDs -> (filename, signature md5)
    signum = 1
    sigid_to_siginfo = {}

    # for every minhash in every signature, link it to its NCBI taxonomic ID.
    if args.traverse_directory:
        inp_files = list(traverse_find_sigs(args.sigs))
    else:
        inp_files = list(args.sigs)

    print('loading signatures & traversing hashes')
    bad_input = 0
    for n, filename in enumerate(inp_files):
        if n % 100 == 0:
            print('... loading file #', n, 'of', len(inp_files), end='\r')

        try:
            sig = sourmash_lib.signature.load_one_signature(filename,
                                                      select_ksize=args.ksize)
        except (FileNotFoundError, ValueError):
            if not args.traverse_directory:
                raise

            bad_input += 1
            continue

        # record which file
        md5sum = sig.md5sum()
        sigid_to_siginfo[signum] = (filename, md5sum)
        this_signum = signum
        signum += 1

        if sig.minhash.scaled < args.scaled:
            sig.minhash = sig.minhash.downsample_scaled(args.scaled)

        # add hashval -> signature info
        mins = sig.minhash.get_mins(with_abundance=True)
        for m, abund in mins.items():
            hashval_to_sigids[m].append(this_signum)
            hashval_to_abunds[m].append(abund)

    print('\n...done. {} hashvals, {} signatures'.format(len(hashval_to_sigids), len(sigid_to_siginfo)))

    if bad_input:
        print('failed to load {} of {} files found'.format(bad_input,
                                                           len(inp_files)))

    hashval_picklefile = args.savename + '.hashvals'
    sig_picklefile = args.savename + '.siginfo'
    abund_picklefile = args.savename + '.abunds'

    with open(hashval_picklefile, 'wb') as hashval_fp:
        dump(hashval_to_sigids, hashval_fp)
    with open(sig_picklefile, 'wb') as sig_fp:
        dump(sigid_to_siginfo, sig_fp)
    with open(abund_picklefile, 'wb') as abund_fp:
        dump(hashval_to_abunds, abund_fp)


if __name__ == '__main__':
    main()
