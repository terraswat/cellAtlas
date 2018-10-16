#!/usr/bin/env bash

LOCALPY=$1

# temporary:
#LOCALPY=/soe/swat/packages/localpy


# To find libs and headers installed:
# yum info zlib zlib-devel

# Build python.
mkdir $LOCALPY
cd $LOCALPY
wget https://www.python.org/ftp/python/3.7.0/Python-3.7.0.tar.xz
tar xf Python-3.7.0.tar.xz
cd Python-3.7.0
./configure --prefix=$LOCALPY
make &> $LOCALPY/Python-3.7.0/make.log
make install &> $LOCALPY/Python-3.7.0/install.log

# Add the bin directory to the path.
PATH=$LOCALPY/bin:$PATH

# Note: pip is built into python3, executed with:
#           python3 -m pip install SomePackage
#       virtual environment is installed by default, executed with:
#           python3 -m venv /path/to/new/virtual/environment

echo
echo 'python3 built in:'
echo $LOCALPY
echo 'python executable:'
which python3



