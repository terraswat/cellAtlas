#!/usr/bin/env bash

export LOCALPY_BIN=/soe/swat/packages/localpy/bin
export SANDBOX=/projects/sysbio/users/cellAtlas
export PYENV=$SANDBOX/env

# Start off pointing to only the newly installed python to minimize path
# problems/confusion
PATH=$LOCALPY_BIN:$PATH

echo LOCALPY_BIN: $LOCALPY_BIN

mkdir $PYENV
python3 -m venv $PYENV

if [ $? != 0 ]; then
   exit $?
fi
# Add bin and source code dirs to the python path in the virtual environment.
# TODO the below path should be more strict by not including PATH from user's
# environment
echo "export PATH=$PYENV/bin:$PATH" >> $PYENV/bin/activate
echo "export LD_LIBRARY_PATH=/soe/swat/packages/hdf5/lib" >> $PYENV/bin/activate
echo "export HDF5_DIR=/soe/swat/packages/hdf5" >> $PYENV/bin/activate
cd $SANDBOX

# Activate the python environment and install the requirements.
source $PYENV/bin/activate # this is not working for some reason
#python3 -m pip install scanpy
pip install -r $SANDBOX/build/requirements.txt


