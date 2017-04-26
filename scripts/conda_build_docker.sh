#!/bin/bash

conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

conda config --set anaconda_upload yes

conda build ./conda
