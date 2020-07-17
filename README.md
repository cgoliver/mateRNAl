# mateRNAl

`mateRNAl` is an RNA evolutionary simulation tool implemented in Python. Simulations can be based on two fitness measurements: base pair distance to a target, or energy of MFE structure on sequence. Simulations can also be constrained to only allow sequences with a specified GC content. 

### Citing

```
@article{oliver2019emergence,
  title={On the emergence of structural complexity in RNA replicators},
  author={Oliver, Carlos G and Reinharz, Vladimir and Waldisp{\"u}hl, J{\'e}r{\^o}me},
  journal={RNA},
  volume={25},
  number={12},
  pages={1579--1591},
  year={2019},
  publisher={Cold Spring Harbor Lab}
}
```
### Requirements

* python 2.7+
* [begins 0.9](https://pypi.python.org/pypi/begins/0.9)
* [Vienna RNA package 2.0+](https://www.tbi.univie.ac.at/RNA/)

### Usage

Print the help menu and options:

```
python mateRNAl.py -h 
```

Default run:

```
python mateRNAl.py 
```
Energy based selection on sequences of length 100 with GC content of 0.3, mutation rate of 0.1, and 2000 generations.

```
python mateRNAl.py -l 100 -g 0.3 -m 0.1 -t 2000
```

When a simulation is to be repeated with the same parameters, the -r flag tells mateRNAl how many replicates to run. Replicate runs can be run in parallel using -p to specify number of processes to use.

```
python mateRNAl.py -l 100 -g 0.3 -m 0.1 -t 2000 -n 5 -p 5
```
### Output

Each run produces a `.csv` file that can be easily parsed with tools such as [pandas](http://pandas.pydata.org/). Each column of the csv is labeled as follows:

```
generation, sequence, structure, energy, probability, gc, mutations, 
```
Manuscript submitted. [Preprint](https://www.biorxiv.org/content/early/2017/11/15/218990)
