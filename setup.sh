#! /bin/bash
#
# file: setup.sh
#
# This bash script performs any setup necessary in order to test your
# entry.  It is run only once, before running any other code belonging
# to your entry.

set -e
set -o pipefail

# Example: compile a C module (viterbi_Schmidt) for Octave
#mkoctfile -mex viterbi_Springer.c

# Example: compile a C module for Matlab
#mex viterbi_Springer.c

# Remove (or set it to 0) if you are not using Matlab
NEED_MATLAB=1

unzip scipy-0.19.0.zip
cd scipy-0.19.0
python setup.py install --user
cd ..

cd Keras-1.2.2
python setup.py install --user
cd ..

# Install matlab-engine for python
MATLAB_BIN=$(dirname $(which matlab))
cd $MATLAB_BIN/../extern/engines/python && python setup.py build --build-base=/tmp/mlpy install --user
