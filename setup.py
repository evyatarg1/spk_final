from setuptools import setup, Extension

setup(name='mykmeanssp',
    version='1.0',
    description='C implementation of kmeans',
    ext_modules=[Extension('mykmeanssp', sources=['spkmeans.c', 'spkmeansmodule.c'])])
