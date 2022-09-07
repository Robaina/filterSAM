# Filtering sam/bam files by percent identity or percent of matched sequence

Tools to filter alignments in SAM/BAM files by percent identity or percent of matched sequence. 

Percent identity is computed as:

<img src="https://render.githubusercontent.com/render/math?math=PI = 100 \frac{N_m}{N_m + N_i}">

where N<sub>m</sub> is the number of matches and N<sub>i</sub> is the number of mismatches.

Percent of matched sequences is computed as:

<img src="https://render.githubusercontent.com/render/math?math=PM = 100 \frac{N_m}{L}">

where L corresponds to query sequence length.

## NOTES

BAM/SAM files must contain [MD tags](https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files) to be able to filter by percent identity. Aligners such as [BWA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2705234/) add MD tags to each queried sequence in a BAM file. MD tags can also be generated with [samtools](http://www.htslib.org/doc/samtools-calmd.html).

## Dependencies

1. [Samtools](http://www.htslib.org/)
2. [Pysam](https://pysam.readthedocs.io/en/latest/api.html)

## Installation

```pip3 install filtersam```

Better to install within an environment, such as a conda environment, to avoid
path conflicts with the included bash scripts.

## TODO

1. Make it command line callable
2. Perhaps good idea (if possible) to add a specific tag to BAM/SAM containing computed percent identity
3. Include several definitions of percent identity and/or let the user define one

# Usage

This package contains two main functions: ```filterSAMbyIdentity``` and ```filterSAMbyPercentMatched```, to filter BAM files by percent identity or percent of matched sequence, respectively. 

To exemplify its usage, let's filter a BAM file by percent identity and percent of matched sequence.


```python
from filtersam.filtersam import filterSAMbyIdentity, filterSAMbyPercentMatched


# Filter alignments with percent identity greater or equal to 95%
filterSAMbyIdentity(input_path='ERS491274.bam',
                    output_path='ERS491274_PI95.bam',
                    identity_cutoff=95)

# Filter alignments with percent of matched sequence greater or equal to 50%
filterSAMbyPercentMatched(input_path='ERS491274.bam',
                          output_path='ERS491274_PM50.bam',
                          matched_cutoff=50)
```

# Parallelizing filtersam

Filtering large BAM files can take a while. However, ```filtersam``` can be parallelized with an additional python package: [parallelbam](https://pypi.org/project/parallelbam/). Effectively, ```parallelbam``` splits a large BAM file into chunks and calls ```filtersam``` in dedicated processes for each one of them.

Let's try this out, we will parallelize the above operation in 8 processes.


```python
from parallelbam.parallelbam import parallelizeBAMoperation, getNumberOfReads


# Filter alignments with percent identity greater or equal to 95% in parallel
parallelizeBAMoperation('ERS491274.bam',
                        callback=filterSAMbyIdentity,
                        callback_additional_args=[95],
                        n_processes=8,
                        output_path='ERS491274_PI95_parallel.bam')
```

We can further check if the filtered bam files produced in a single process and in parallel contain the same number of segments with the function ```getNumberOfReads``` of parallelbam.


```python
# Number of segments in the original bam
getNumberOfReads('ERS491274.bam')
```




    1113119




```python
# Number of segments in the single-process PI-filtered bam file
getNumberOfReads('ERS491274_PI95.bam')
```




    11384




```python
# Number of segments in the paralllized PI-filtered bam file
getNumberOfReads('ERS491274_PI95_parallel.bam')
```




    11384



We see that both bam files contain the same number of (filtered) segments (fewer than in the original bam file).

# Command-line usage

Filtersam can also be called as a command line program in the following way:

```filtersam [-h] [-i] [-m] [-p] [-o] bam```

where _bam_ is the path to the bam/sam file.

Call 

```filtersam --help```

to display help text about the arguments.
