'''hmMHC - a hidden Markov model-based MHC binding predictor

Copyright (c) 2019 Maxim Artyomov, Ilya Kizhvatov
'''

from __future__ import division
import numpy as np
import pandas as pd
import ghmm
import logging
import pkg_resources

# suppress excessive WARNING-level messages from ghmm.loglikelihoods()
logging.getLogger("GHMM").setLevel(logging.ERROR)

class hmMHC:
    ''' A hidden Markov model-based MHC binding predictor '''

    def __init__(self, allele):
        ''' Initialize the predictor. For now, only H2-IAb is supported. '''

        if (allele != 'H2-IAb'):
            raise ValueError('Unsupported allele')

        self.allele = allele

        # range of allowed peptide lengths (percentile rank calibration was
        # done for this range)
        self.minPeptideLength = 12
        self.maxPeptideLength = 24

        # define the HMM emission alphabet, including a dedicated "stop"
        # emission (a GHMM-specific tweak)
        self.alphabet = ghmm.Alphabet(ghmm.AminoAcids.listOfCharacters + ['Z'])

        # a value to replace inf in the predictions. Same value as used in
        # percentile rank calibration
        self.infSubstitute = 1000

        # GHMM has an internal limitation of 1,500,000 sequences in a set
        self.sequencesPerBlock = 1500000

        # load the HMM model and the corresponding percentile rank calibration data
        from pkg_resources import resource_filename
        hmmFile = resource_filename(__name__, 'models/hmm-h2iab-iedb2018-binders.xml')
        prCalibrationFile = resource_filename(__name__, 'models/precent-rank-model-h2iab.npz')

        self.hmm = ghmm.HMMOpen(hmmFile)

        prCalibration = np.load(prCalibrationFile)
        self.prCalibrationCdf = prCalibration['cdf']
        self.prCalibrationBinEdges = prCalibration['bin_edges']

    def predict(self, peptides):
        ''' Computes binding predictions for a list of peptides.

        Returns a dataframe with peptides, predicted -10 log odds
        values, and calibrated percentile ranks
        '''

        # accept both list of strings and a single peptide as a string
        if isinstance(peptides, str):
            peptidesList = [peptides]
        elif isinstance(peptides, list):
            peptidesList = peptides
        else:
            raise ValueError('Wrong input format, must be a list or a single string')

        # predict -10 length-normalized log odds values
        normalizedLogOdds = self.computeLogOdds(peptidesList)

        # convert -10 log odds to calibrated percentile rank
        percentileRanks = self.computePercentileRanks(normalizedLogOdds)

        # combine everything into the table (the order is preserved)
        return pd.DataFrame(
            zip(peptidesList, normalizedLogOdds, percentileRanks),
            columns = ['peptide', '-10logOdds', 'percentile_rank']
        )

    def computeLogOdds(self, peptides):
        ''' Returns length-normalized log-odds values from the model for every 
        peptide in the list, preserving the order
        '''

        # get peptide lengths and check that they are within the working range
        peptideLengths = np.array([len(p) for p in peptides])
        if (
            (peptideLengths.min() < self.minPeptideLength) |
            (peptideLengths.max() > self.maxPeptideLength)
            ):
            raise ValueError(
                'Peptide length outside of the working range {}-{}'.format(
                        self.minPeptideLength,
                        self.maxPeptideLength
                    )
                )

        # ensure we have a list of non-unicode strings and append the "stop" emission
        peptideList = [str(p) + 'Z' for p in peptides]

        # convert to blocks of sequence sets
        sequenceSetBlocks = self.toSequenceSetBlocks(peptideList)

        # get log likelihoods from HMM, per block
        logLikelihoodBlocks = []
        for s in sequenceSetBlocks:
            logL = self.hmm.loglikelihoods(s)
            logLikelihoodBlocks.append(np.array(logL))

        # concatenate the blocks
        logLikelihoods = np.concatenate(logLikelihoodBlocks)
        
        # normalize predictions by peptide length (including the "stop" emission)
        # and scale by -10 such that smaller values correspond to higher binding
        # likelihood (same as with affinity values)
        normalizedLogOdds = -10 * logLikelihoods / (peptideLengths + 1)
        
        # replace inf's by very large values
        normalizedLogOdds[np.isinf(normalizedLogOdds)] = self.infSubstitute

        return normalizedLogOdds

    def computePercentileRanks(self, normalizedLogOdds):
        ''' Compute percentile ranks for a list of length-normalized log odds values '''
        
        indices = np.searchsorted(self.prCalibrationBinEdges, normalizedLogOdds)
        percentileRanks = self.prCalibrationCdf[indices]

        return percentileRanks

    def toSequenceSetBlocks(self, peptideList):
        ''' Converts a list of peptides given as strings to GHMM format.

        As GHMM has a limitation of max 1,500,000 sequences per sequence set,
        longer peptide lists are split into blocks. The block size is
        configurable via self.sequencesPerBlock to facilitate testing.

        Returns a list of sequence sets, preserving the original order
        of the peptides.
        '''
        
        # split into blocks 
        lenPeptides = len(peptideList)
        numFullBlocks = lenPeptides // self.sequencesPerBlock
        lenLastBlock = lenPeptides % self.sequencesPerBlock
        
        sequenceSets = []

        # full blocks (if any)
        for i in range(0, numFullBlocks):
            rangeStart = i * self.sequencesPerBlock
            rangeEnd = rangeStart + self.sequencesPerBlock
            sequenceBlock = ghmm.SequenceSet(
                self.alphabet,
                [list(p) for p in peptideList[rangeStart:rangeEnd]]
            )
            sequenceSets.append(sequenceBlock)

        # the last partial block (if any)
        rangeStart = numFullBlocks * self.sequencesPerBlock
        rangeEnd = lenPeptides
        if (lenLastBlock > 0):
            sequenceBlock = ghmm.SequenceSet(
                    self.alphabet,
                    [list(p) for p in peptideList[rangeStart:rangeEnd]]
                )
            sequenceSets.append(sequenceBlock)
        
        return sequenceSets
