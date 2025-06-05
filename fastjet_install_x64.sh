#!/bin/bash

setup 64b
FASTJET_VERSION=3.4.2
FASTJET_CONTRIB_VERSION=1.046

install_dir=/gpfs01/star/pwg/$USER/install/x64
current_dir=$(pwd)

# install fastjet
mkdir -p $install_dir
cd $install_dir
curl -O https://fastjet.fr/repo/fastjet-$FASTJET_VERSION.tar.gz
tar zxvf fastjet-$FASTJET_VERSION.tar.gz
cd fastjet-$FASTJET_VERSION/

./configure --prefix=$PWD/../fastjet-install
make
make check
make install

#  install fastjet-contribcd
cd ..
wget http://fastjet.hepforge.org/contrib/downloads/fjcontrib-$FASTJET_CONTRIB_VERSION.tar.gz
tar zxvf fjcontrib-$FASTJET_CONTRIB_VERSION.tar.gz
cd fjcontrib-$FASTJET_CONTRIB_VERSION
./configure --fastjet-config=$PWD/../fastjet-$FASTJET_VERSION/fastjet-config --prefix=$PWD/../fastjet-install
make
make install

# add FASTJET and FASTJET_CONTRIB to /.bashrc
export FASTJET=$install_dir/fastjet-install
export FASTJET_CONTRIB=$install_dir/fjcontrib-$FASTJET_CONTRIB_VERSION

echo "export FASTJET='$install_dir/fastjet-install'" >>~/.bashrc
echo "export FASTJET_CONTRIB='$install_dir/fjcontrib-$FASTJET_CONTRIB_VERSION'" >>~/.bashrc

source ~/.bashrc
