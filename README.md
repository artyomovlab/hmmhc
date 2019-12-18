# hmMHC 
A hidden Markov model-based MHC II binding predictor. Currently, only H2-IAb predictions are supported. The predictor is described in the paper

Elise Alspach et al. MHC-II neoantigens shape tumour immunity and response to immunotherapy. Nature (2019), available at https://www.nature.com/articles/s41586-019-1671-8

## Installation
hmMHC can be installed on Linux and macOS via a combination of conda and pip install. Windows is not supported. Python 3 is not supported.
```
$ conda create -n hmmhc -c bioconda python=2.7 ghmm=0.9 'icu=58.*'
$ conda activate hmmhc
$ pip install git+https://github.com/artyomovlab/hmmhc#egg=hmmhc
```

## Command line example
Input from command line and output to stout:
```
$ hmmhc-predict --allele H2-IAb --peptides VNGYNEAIVHVVETP IKSEHPGLSIGDVAK KESVVSGKAVPREEL
```
Input from csv and output to csv:
```
$ hmmhc-predict --allele H2-IAb --input example_input.csv --output output.csv
```
See `hmmhc-predict -h` for further details.

## Python example
```python
from hmmhc import hmMHC
predictor = hmMHC('H2-IAb')

peptides = ['VNGYNEAIVHVVETP', 'IKSEHPGLSIGDVAK', 'KESVVSGKAVPREEL']

predictor.predict(peptides)
```

## Output

The predictor outputs a list of peptides with the predicted -10 log odds scores and corresponding percentile ranks. Percentile ranks are computed from -10 log odds scores based on model calibration on a large set of random natural peptides. For both metrics, smaller values correspond to higher binding likelihood. See Methods section in the paper for further details.

## Dependencies

hmMHC relies on General Hidden Markov Model library (GHMM) by A. Schliep et al., see http://ghmm.sourceforge.net/.

## Latest version

The latest version of hmMHC is available at https://github.com/artyomovlab/hmmhc.
