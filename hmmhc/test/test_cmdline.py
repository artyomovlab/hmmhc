'''Tests for hmMHC command line interface'''

from hmmhc.cmdline import main
import io
import sys
import tempfile
import pandas as pd
from nose.tools import *


# universal test vector
inList = ['VNGYNEAIVHVVETP', 'IKSEHPGLSIGDVAK', 'KESVVSGKAVPREEL']
inCsv = '''
VNGYNEAIVHVVETP
IKSEHPGLSIGDVAK
KESVVSGKAVPREEL
'''
expectedCsvOutput = '''
peptide,-10logOdds,percentile_rank
VNGYNEAIVHVVETP,26.096784479388965,1.7982905702723602
IKSEHPGLSIGDVAK,28.51906056754387,14.708167060725673
KESVVSGKAVPREEL,26.195585528673668,1.9979640447254532
'''

# command line input and output to stdout
@raises(SystemExit)
def testCommandLine():

    sysStdout = sys.stdout
    out = io.BytesIO()
    sys.stdout = out
    
    main(['--allele', 'H2-IAb', '--peptides'] + inList)
    
    assert out == expectedCsvOutput
    sys.stdout = sysStdout

# command line input and output to stdout
@raises(SystemExit)
def testCommandLineCsv():

    inFile = tempfile.NamedTemporaryFile(suffix='.csv')
    inFile.write(inCsv)
    inFile.flush()
    outFile = tempfile.NamedTemporaryFile(suffix='.csv')
    
    main(['--allele', 'H2-IAb', '--input', inFile.name, '--output', outFile.name])

    res = outFile.read()

    assert res == expectedCsvOutput

    inFile.close()
    outFile.close()
