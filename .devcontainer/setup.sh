#!/bin/bash

# OS packages
## apt-get update
apt-get update -y
## apt-get install
apt-get install -y build-essential tree screen git emacs r-base

# conda env install
micromamba env create -y -f conda_env.yaml -n phyloflash
