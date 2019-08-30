# hmMHC 
A hidden Markov model-based MHC II binding predictor. Currently, only H2-IAb predictions are supported.

## Installation
hmMHC can be installed on Linux and macOS via a combination of conda and pip install. Windows is not supported. Python 3 is not supported.
```
$ conda create -n hmmhc -c bioconda python=2.7 ghmm=0.9 'icu=58.*'
$ conda activate hmmhc
$ pip install git+https://github.com/artyomovlab/hmmhc#egg=hmmhc
```

## Command line example
```
$ hmmhc-predict --allele H2-IAb --peptides VNGYNEAIVHVVETP IKSEHPGLSIGDVAK KESVVSGKAVPREEL --out predictions.csv
```
See `hmmhc-predict -h` for further details.

## Python example
```python
from hmmhc import hmMHC
predictor = hmMHC('H2-IAb')

peptides = ['VNGYNEAIVHVVETP', 'IKSEHPGLSIGDVAK', 'KESVVSGKAVPREEL']

predictor.predict(peptides)
```

## Dependencies

hmMHC relies on General Hidden Markov Model library (GHMM) by A. Schliep et al., see http://ghmm.sourceforge.net/.
