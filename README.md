blockTools
=========

*A bSFS-based analysis toolkit*

Dependencies (via [conda](https://conda.io/miniconda.html))
-------

```
# clone repository
git clone https://github.com/DRL/blocktools.git

# create conda enviroment with dependencies
conda env install -f $BLOCKTOOLS_PATH/blocktools.conda.yml

# Activate blocktools conda environment
conda activate blocktools
```

Usage
-----

- Input files are specified in the ```yaml``` files (see ```data/input/toy.multibed.A.yaml``` and ```data/input/toy.multibed.B.yaml```)
- Algorithm A is faster than B, but algorithm B makes more (divergent) blocks

```
# Making blocks
./blocktools blocks -y data/input/toy.multibed.A.yaml

# Fetching variants from VCF
./blocktools variants -y data/input/toy.multibed.A.yaml

# Making windows of blocks (window size = 3 blocks, window overlap = 1 block; real datasets would use -w 500; -l 100)
./blocktools windows -y data/input/toy.multibed.A.yaml -w 3 -l 1

# Making plots
./blocktools plots -y data/input/toy.multibed.A.yaml -g data/input/toy/toy.genomefile
```
