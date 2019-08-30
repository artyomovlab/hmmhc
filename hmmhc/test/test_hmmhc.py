'''Tests for the hmMHC class'''

from hmmhc import hmMHC

from numpy import testing
from nose.tools import *

# A selection of screened mutant peptides from the Nature paper
# as a universal test vector
#
# protein   peptide sequence    -10logOdds  percentile_rank
# Itgb1     VNGYNEAIVHVVETP     26.0967845   1.7982906
# Hmgb1     IKSEHPGLSIGDVAK     28.5190606  14.7081671
# Akr1a1    KESVVSGKAVPREEL     26.1955855   1.9979640
# Zzz3      VDHDADFQGAKPACR     27.6227440   7.4267293
# Eif4a3    KEFRSAASRVLISTD     26.4548705   2.5557151

peptides = [
    'VNGYNEAIVHVVETP', 'IKSEHPGLSIGDVAK', 'KESVVSGKAVPREEL',
    'VDHDADFQGAKPACR', 'KEFRSAASRVLISTD'
]
logodds = [26.0967845, 28.5190606, 26.1955855, 27.6227440, 26.4548705]
pranks = [1.7982906, 14.7081671, 1.9979640, 7.4267293, 2.5557151]


# prediction produces the expected values in the expected order
def testPrediction():

    # run predictor
    predictor = hmMHC('H2-IAb')
    res = predictor.predict(peptides)

    # check
    assert peptides == res.peptide.to_list()
    testing.assert_almost_equal(res['-10logOdds'], logodds)
    testing.assert_almost_equal(res['percentile_rank'], pranks)


# a single peptide is accepted as a string
def testSingleStringInput():

    p = 'VNGYNEAIVHVVETP'

    res = hmMHC('H2-IAb').predict(p)

    assert res.loc[0, 'peptide'] == p


#  peptide legnth limit is correctly enforced
def testAdmissibleLengths():

    # peptides with all admissible lengths, from 12 to 24
    peptidesValid = [
        'GFAVVRPPGHHA', # 12
        'LISTARPSFMDLP',
        'NRQPPSVRPNQHHF',
        'AERSAAQESAHLGGP',
        'GISPGDSTTNDAPHSG',
        'PPEDRNSVAAMQSEPGS',
        'SKGNIASQKSDYLKHCTF',
        'QENIKVLESDLSEEREKRQ',
        'PTITPLVQISSDKRIINVLK',
        'SPYFYPQSLVSNLDPGAALYL',
        'DCHRQLKDSKQILSITKNFKVE',
        'GKQHSESLNIRVYQPPAQVTLKL',
        'DIVEVLFTQPNVELNQQNKLGDTA' # 24
    ]

    predictor = hmMHC('H2-IAb')

    # this should execute without exceptions
    predictor.predict(peptidesValid)


# peptide lower length limit is correctly enforced
@raises(ValueError)
def testLowerLengthLimt():

    # a short peptide is present
    peptidesShort = [
        'VVKGIRLSENVI', # 12
        'KEKKWSKVGSR', # 11
        'TDMKDKRNLTEFRQLVYCSAVKNF' # 24
    ]

    predictor = hmMHC('H2-IAb')

    # this should raise a ValueError exception
    predictor.predict(peptidesShort)


# peptide upper length limit is correctly enforced
@raises(ValueError)
def testUpperLengthLimt():

    # a long peptide is present
    peptidesLong = [
        'EANGSTAWPPPTASNISEPHQCLL', # 24
        'TVPKAGTVPLATEVLKNLTAPPTLE', # 25
        'WMLMAELGTIET' # 12
    ]

    predictor = hmMHC('H2-IAb')

    # this should raise a ValueError exception
    predictor.predict(peptidesLong)


# block splitting works as expected, preserving the order
# same as testPrediction, but block size is reduced to trigger splitting
def testBlockSplitting():

    predictor = hmMHC('H2-IAb')
    predictor.sequencesPerBlock = 2
    res = predictor.predict(peptides)

    assert peptides == res.peptide.to_list()
    testing.assert_almost_equal(res['-10logOdds'], logodds)
    testing.assert_almost_equal(res['percentile_rank'], pranks)


# block splitting is actually performed
def testBlockSplittingLowLevel():

    predictor = hmMHC('H2-IAb')

    # split a list of 5 peptides into blocks of max len 2
    predictor.sequencesPerBlock = 2
    listOfBlocks = predictor.toSequenceSetBlocks(peptides)

    # check that there are 3 blocks of lengths 2, 2, 1
    assert len(listOfBlocks) == 3
    assert len(listOfBlocks[0]) == 2
    assert len(listOfBlocks[1]) == 2
    assert len(listOfBlocks[2]) == 1

    # split a list of 4 peptides into blocks of max len 2
    predictor.sequencesPerBlock = 2
    listOfBlocks = predictor.toSequenceSetBlocks(peptides[0:4])

    # check that there are 2 blocks of lengths 2, 2
    assert len(listOfBlocks) == 2
    assert len(listOfBlocks[0]) == 2
    assert len(listOfBlocks[1]) == 2
