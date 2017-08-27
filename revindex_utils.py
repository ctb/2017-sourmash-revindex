#! /usr/bin/env python
from pickle import load

class HashvalRevindex(object):
    def __init__(self, savename):
        x = load_revindex(savename)
        self.hashval_to_sigids, self.sigid_to_siginfo, self.hashval_to_abunds = x


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