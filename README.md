![logo](assets/logo.png)
## A Python tool to filter sam/bam files by percent identity or percent of matched sequence

[![DOI](https://zenodo.org/badge/400865776.svg)](https://zenodo.org/badge/latestdoi/400865776)

<br>

Percent identity is computed as:

$$PI = 100 \frac{N_m}{N_m + N_i}$$

where $N_m$ is the number of matches and $N_i$ is the number of mismatches.

Percent of matched sequences is computed as:

$$PM = 100 \frac{N_m}{L}$$

where $L$ corresponds to query sequence length.

## NOTE

BAM/SAM files must contain [MD tags](https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files) to be able to filter by percent identity. Aligners such as [BWA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2705234/) add MD tags to each queried sequence in a BAM file. MD tags can also be generated with [samtools](http://www.htslib.org/doc/samtools-calmd.html).

## Installation

```pip install filtersam```

## Usage

You can find a jupyter notebook with usage examples [here](examples/examples.ipynb).

## Citation

If you use this software, please cite it as below:

Robaina-Estévez, S. (2022). filterSAM: filter sam/bam files by percent identity or percent of matched sequence (Version 0.0.11)[Computer software]. https://doi.org/10.5281/zenodo.7056278.