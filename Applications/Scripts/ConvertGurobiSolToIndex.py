#!/usr/bin/env python
"""
Description: Covert .sol file (by gurobi) to index file for subsampling gradients.
It is used for subsampling gradients by Mixed Integer Linear Programming (MILP)

Usage:
  ConvertGurobiSolToIndex.py  <solfile> <indexfile> -n <numberOfSamples> [-v]
  ConvertGurobiSolToIndex.py (-h | --help)

Options:
  -h --help                Show this screen.
  -v --verbose             Verbose
  -n <numberOfSamples>     number of samples used in shells

Examples:
ConvertGurobiSolToIndex.py solution_gurobi.sol indexInShells.txt -n 90,90,90

Author(s): Jian Cheng (jian.cheng.1983@gmail.com)
"""

import os
from docopt import docopt


def main():

    args = docopt(__doc__, version='1.0')

    if (args['--verbose']):
        print(args)

    _solfile = args['<solfile>']
    _indexfile = args['<indexfile>']
    _verbose = args['--verbose']
    _num = args['-n']

    with open(_solfile) as f:
        sol = f.readlines()

    sol = [x for x in sol if x.strip() and x[0]!='#']
    if _verbose:
        print sol

    numbers = [int(x) for x in _num.split(',')]
    numerOfShells = len(numbers)
    if sum(numbers)+numerOfShells+1 != len(sol):
        raise Exception("wrong number of samples!\n")

    index_main, index_ext = os.path.splitext(_indexfile)

    numberDoubles = numerOfShells+1 if numerOfShells>1 else numerOfShells
    for s in xrange(0,numerOfShells):
        i_start = numberDoubles + sum(numbers[:s])
        #  print s, i_start
        index_list = []
        for i in xrange(i_start, i_start+numbers[s]):
            binary =  int(sol[i].split(' ')[1].strip())
            if binary==1:
                index_list.append(i-i_start)
        #  print index_list, len(index_list)

        file_name = _indexfile if numerOfShells==1 else ('%s_shell%d%s' % (index_main, s+1, index_ext))
        print 'write index file of shell %d to %s' % (s+1, file_name)
        with open(file_name, mode='wt') as f:
            f.write('\n'.join(str(ind) for ind in index_list))


if __name__ == '__main__':
    main()
