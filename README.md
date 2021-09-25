TAU Software Project Course - Final project

# Spectal-Clustering
Implementation of K-means and Spectral Clustering algorithm

Submitted by:
- Evyatar Gavish
- Daniella Felig

Compile on C:
>>>gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans
Run in C:
>>>./spkmeans K Goal input.txt 

Python setup:
>>>python3 setup.py build_ext --inplace
Run in Python:
>>>python3 spkmeans.py 3 spk input.txt
*k is int smaller then n; goal is "wam", "ddg", "lnorm", "jacobi" of "spk".

Relevant files:
1. spkmeans.py: Python interface of your code.
2. spkmeans.h: C header file.
3. spkmeans.c: C interface of your code.
4. spkmeansmodule.c: Python C API wrapper.
5. setup.py: The setup file.
