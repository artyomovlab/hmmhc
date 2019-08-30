'''hmMHC command line interface

Copyright (c) 2019 Maxim Artyomov, Ilya Kizhvatov
'''

from __future__ import print_function
import argparse
import pandas as pd
import sys
from .hmmhc import hmMHC


def parseArgs(args): 
    '''Define and parse command line arguments'''

    parser = argparse.ArgumentParser(description='hmMHC - a hidden Markov model-based MHC binding predictor')

    parser.add_argument('--allele', help='allele (currently, only H2-IAb is supported)', type=str, required=True)
    parser.add_argument('--output', help='output CSV file name', type=str, metavar='FILENAME')

    inputArgGroup = parser.add_mutually_exclusive_group()
    inputArgGroup.add_argument('--input', help='input CSV file name (exclusive with --peptides)', type=str, metavar='FILENAME')
    inputArgGroup.add_argument(
        '--peptides',
        help='peptide sequences delimited by whitespaces (exclusive with --in)',
        nargs = '+',
        type=str,
        metavar='PEPTIDE'
    )
    
    return parser.parse_args(args)


def main(arguments=sys.argv[1:]):
    '''hmMHC command line entrypoint'''

    # get command line arguments
    args = parseArgs(arguments)

    # get peptides from input
    if args.input:
        peptidesDf = pd.read_csv(args.input, header=None)
        peptides = peptidesDf[0].to_list()
    elif args.peptides:
        peptides = args.peptides
    else:
        parser.print_usage(sys.stderr)
        print('Error: no input provided', file=sys.stderr)
        exit(1)

    # predict
    predictor = hmMHC('H2-IAb')
    predictions = predictor.predict(peptides)

    # output predictions
    if (args.output):
        predictions.to_csv(args.output, index=False)
    else:
        predictions.to_csv(sys.stdout, index=False)

    exit(0)
