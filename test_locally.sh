#!/bin/bash

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh

for TRAVIS_PYTHON_VERSION in 3.6 3.7 3.8; do
    conda create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda -c conda-forge python=$TRAVIS_PYTHON_VERSION --file requirements_travis.txt
    conda activate pgt_test_${TRAVIS_PYTHON_VERSION}
    conda install --yes -c conda-forge -c bioconda pytest ghostscript coverage coverage-badge
    python setup.py install
    coverage run -m py.test
    coverage html
    coverage-badge -f -o docs/coverage.svg
    conda deactivate
done

TRAVIS_PYTHON_VERSION=3.9
conda create -n pgt_test_${TRAVIS_PYTHON_VERSION} --yes -c bioconda -c conda-forge python=$TRAVIS_PYTHON_VERSION
conda activate pgt_test_${TRAVIS_PYTHON_VERSION}
conda install --yes -c conda-forge -c bioconda numpy bedtools h5py pytest ghostscript
pip install -r requirements_travis.txt
python setup.py install
py.test pygenometracks --doctest-modules
conda deactivate
