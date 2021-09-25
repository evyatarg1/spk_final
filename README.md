# TAU Software Project Course - Final project

Spectal-Clustering
Implementation of K-means and Spectral Clustering algorithm

Submitted by:
- Evyatar Gavish
- Daniella Felig

Compile on C:
>>>gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans

Run in C:
>>>./spkmeans K Goal input.txt 

Reading user CMD arguments:

(a) k (int, < N): Number of required clusters.

(b) goal : Can get the following values:
- spk: results of the normalized spectral clustering algorithm
- wam: Weighted Adjacency Matrix
- ddg: Diagonal Degree Matrix
- lnorm: Normalized Graph Laplacian
- jacobi: eigenvalues and eigenvectors of the points being represented as a real symmetric matrix 
  
(c) file name (.txt or .csv): The path to the file that will contain N observations

Python setup:
>>>python3 setup.py build_ext --inplace

Run in Python:
>>>python3 spkmeans.py 3 spk input.txt

(Same instructions as C part)

Relevant files:
1. spkmeans.py: Python interface of your code.
2. spkmeans.h: C header file.
3. spkmeans.c: C interface of your code.
4. spkmeansmodule.c: Python C API wrapper.
5. setup.py: The setup file.
