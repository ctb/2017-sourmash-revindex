#! /usr/bin/env python
# @CTB we should be storing the ksize here, although we can get it from
# the signatures...
from pickle import load
import sourmash_lib.signature

class HashvalRevindex(object):
    def __init__(self, savename):
        x = load_revindex(savename)
        self.hashval_to_sigids, self.sigid_to_siginfo, self.hashval_to_abunds = x

    def get_siginfo(self, hashval):
        for siginfo_id in self.hashval_to_sigids.get(hashval, []):
            yield self.sigid_to_siginfo[siginfo_id]


def get_sourmash_signature(filename, md5):
    "Get sourmash signature based on siginfo."
    sigiter = sourmash_lib.signature.load_signatures(filename)
    for sig in sigiter:
        if sig.md5sum() == md5:
            return sig

    raise Exception('no such sig: {}:{}'.format(filename, md5))


def load_revindex(savename):
    hashval_picklefile = savename + '.hashvals'
    sig_picklefile = savename + '.siginfo'
    abunds_picklefile = savename + '.abunds'

    with open(hashval_picklefile, 'rb') as hashval_fp:
        hashval_to_sigids = load(hashval_fp)
    with open(sig_picklefile, 'rb') as sig_fp:
        sigid_to_siginfo = load(sig_fp)
    with open(abunds_picklefile, 'rb') as abund_fp:
        hashval_to_abunds = load(abund_fp)

    return (hashval_to_sigids, sigid_to_siginfo, hashval_to_abunds)


if __name__ == '__main__':
    import sys
    print('loading', sys.argv[1])
    h = HashvalRevindex(sys.argv[1])
