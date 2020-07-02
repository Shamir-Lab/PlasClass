# PlasClass
This module allows for easy classification of sequences as either plasmid or chromosomal.
For example, it can be used to classify the contigs in a (metagenomic) assembly.

## Installation

`plasclass` is written in Python3 and requires NumPy and scikit-learn and their dependencies. These will be installed by the setup.py script.

We recommend using a virtual environment. For example, in Linux, before running setup.py:
```
python -m venv classification-env
source classification-env/bin/activate
```
In Windows:
```
pip install virtualenv
virtualenv classification-env
classification-env\Scripts\activate
```

To install, download and run setup.py:

    git clone https://github.com/Shamir-Lab/PlasClass.git
    cd PlasClass
    python setup.py install

It is possible to install as a user without root permissions:
```
python setup.py install --user
```

After installing, run the tests:
```
python test/test.py
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
my_classifier = plasclass.plasclass()
my_classifier.classify(seqs)
```
The `plasclass()` constructor takes optional parameters:

`n_procs` - number of processes to use for classification. Default=1.

`scales` - array of the scales for the sequence lengths. Default=[1000,10000,100000,500000]

`ks` - array of the k-mer lengths. Default=[3,4,5,6,7]

The sequence(s) to classify, `seqs`, can be either a single string or a list of strings. The strings must be uppercase.

The function `plasclass.classify(seqs)` returns a list of plasmid scores, one per input sequence, in the same order as the input.

### Training new models

The script `train.py` can be used to train new models:
```
python train.py -p <plasmid file> -c <chromosome file> -o <output directory> [-n <num processes> default: 16] [-k <kmer lengths> default: 3,4,5,6,7] [-l <sequence lengths> default: 1000,10000,100000,500000]
```
The command line options for this script are:

`-p/--plasmid`: The fasta file of the plasmid references.

`-c/--chromosome`: The fasta file of the chromosome references.

`-n/--num_processes`: Number of processes to use.

`-o/--outdir`: The path of the output directory. Default=`bin`.

`-k/--kmers`: Comma separated list of the k-mer sizes to use. Default=3,4,5,6,7.

`-l/--lengths`: Comma separated list of the sequence lengths to use. Default=1000,10000,100000,500000.

The models should be put into the `data` directory.

Note that if k-mer and sequence lengths other than the default are used, then these must be specified when calling the `plasclass()` constructor.
