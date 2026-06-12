# Algorithms in Computational Biology

Implementations of computational biology algorithms built from scratch in Python.

## Algorithms

### Exact Pattern Matching
**`z_algo.py`** — Z-Algorithm for DNA pattern matching

Finds all occurrences of a pattern p in text t using the Z-algorithm. Reports match positions and character comparison statistics.

### Pairwise Sequence Alignment
**`global_align.py`** — Needleman-Wunsch Global Alignment

Global alignment of protein sequences using the BLOSUM62 scoring matrix with a fixed indel penalty of σ = 5.

**`local_align.py`** — Smith-Waterman Local Alignment with Affine Gap Penalty

Local alignment of protein sequences using the BLOSUM62 scoring matrix with affine gap penalty (opening = 11, extension = 1).

## Tech Stack
- Python
- BLOSUM62 scoring matrix
- Dynamic programming

## How to Run

```bash
# Z-Algorithm
python z_algo.py input_file.txt

# Global Alignment
python global_align.py input_file.txt

# Local Alignment
python local_align.py input_file.txt
```

## Author
Aleah Hassabo
