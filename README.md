# PlasClass
This module allows for easy classification of sequences as either plasmid or chromosomal.
For example, it can be used to classify the contigs in a (metagenomic) assembly.

## Installation

`classification` is written in Python2.7 and requires NumPy and scikit-learn (versions compatible with Python2.7) and their dependencies. These will be installed by the setup.py script.

We recommend using a virtual environment. For example, before running setup.py:
```
python -m virtualenv classification-env
source classification-env/bin/activate
```

To install, download and run setup.py:

    git clone https://github.com/Shamir-Lab/PlasClass.git
    cd PlasClass
    python setup.py install

If not using a virtual environment, it is possible to install as a user without root permissions:
```
python setup.py install --user
```


<!--- `classification` can also be installed using `pip`. Just do `pip install classification` --->


## Usage

The script `classify_fasta.py` can be used to classify the sequences in a fasta file:
```
python classify_fasta.py -f <fasta file> [-o <output file> default: <fasta file>.probs.out] [-p <num processes> default: 8]
```
The command line options for this script are:

`-f/--fasta`: The fasta file to be classified

`-o/--outfile`: The name of the output file. If not specified, \<input filename\>.probs.out

`-p/--num_processes`: The number of processes to use. Default=8

The output file is a tab separated file with each line containing a sequence header and the corresponding score. The sequences are in the same order as in the input fasta file. 

The classifier can also be imported and used directly in your own python code. For example, once the `plasclass` module has been installed you can use the following lines in your own code:
```
from plasclass import plasclass
my_classifier = plasclass()
my_classifier.classify(seqs)
```
The `plasclass()` constructor takes an optional parameter of the number of processes to use for classification. The default is 1.

The sequence(s) to classify, `seqs`, can be either a single string or a list of strings. The strings must be uppercase.

The function `plasclass.classify(seqs)` returns a list of plasmid scores, one per input sequence, in the same order as the input.

