# Parallel PTAS for Finding Compact Structural Motifs
This repository contains the implementation of a parallel polynomial-time approximation scheme for finding Compact Structural Motifs (as defined by Qian, et al. 2007[^1]).

## Abstract
Structural motifs refer to patterns in 3D space that bear biological significance as they can indicate regions in the protein that have important roles in biochemical functions. Structural motif-finding methods often suffer from a trade off between speed and accuracy, with better quality solutions often taking a significant amount of time to process. This paper utilizes parallelization techniques on the polynomial-time approximation scheme (PTAS)-based (R,C)-compact structural motif finding algorithm in order to minimize this trade off, empirical results presenting a speedup between 4x - 5x from sequential to parallel using three protein data sets. The results of this study emphasize the potential of such solutions in order to maximize both time and accuracy, especially in the case of PTAS problems.

## How to use
1. Open the file `parallel_rccmp.py`.
2. To choose a dataset folder, go to line 216 and change the path to: `structs = get_structures("/datasets","/<your dataset folder here>/")`
3. Variables for length of motif, sample size, and max. ball size (`BENCHMARK_LENGTH`, `r`, and `b`, respectively) can be found and edited on lines 211 - 213.
4. Executing the file should give you the best/minimal sample, RMSD, occurence, and ball size (`SAMPLE`, `BEST_RMSD`, `BEST_OCCURENCE`, and `BALL SIZE`, respectively), as well as the amount of time the program took in seconds.


## Dependencies
This implementation uses Python ver. 3.9.13.
Necessary downloadable libraries are the following:
- Biopython ver. 1.79[^2]
- tqdm ver. 4.64[^3]
- numpy[^4]

### Authors
This study was authored by Bernard Brocka and Sharlene Yap in fulfilment of their CS199 course under the guidance of their adviser, Jhoirene B. Clemente, Ph.D. 

### References
[^1]: Jianbo Qian, Shuaicheng Li, Dongbo Bu, Ming Li, and Jinbo Xu. Finding compact structural motifs. In Combinatorial Pattern Matching, pages 142–149, 07 2007.
[^2]: Biopython. http://biopython.org/, 2021.
[^3]: tqdm. https://github.com/tqdm/tqdm, 2022.
[^4]: Numpy. https://numpy.org/, 2022.
