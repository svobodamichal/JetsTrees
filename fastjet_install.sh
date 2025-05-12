#!/bin/bash
install_dir=/gpfs01/star/pwg/$USER/install
current_dir=$(pwd)

# install fastjet
mkdir -p $install_dir
cd  $install_dir
curl -O https://fastjet.fr/repo/fastjet-3.4.2.tar.gz 
tar zxvf fastjet-3.4.2.tar.gz
cd fastjet-3.4.2/
./configure --prefix=$PWD/../fastjet-install CXXFLAGS="-m32 -fPIC -fno-inline" CFLAGS="-m32 -fPIC -fno-inline" LDFLAGS="-m32"
make
make check
make install

#  install fastjet-contrib
cd  $install_dir
wget http://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.046.tar.gz
tar zxvf fjcontrib-1.046.tar.gz
cd fjcontrib-1.046
./configure --fastjet-config=$PWD/../fastjet-3.4.2/fastjet-config --prefix=$PWD/../fastjet-install CXXFLAGS="-m32 -fPIC -fno-inline" CFLAGS="-m32 -fPIC -fno-inline" LDFLAGS="-m32"
make
make install
make fragile-shared
make fragile-shared-install
# unset not needed flags
unset CXXFLAGS CFLAGS LDFLAGS

# add FASTJET and FASTJET_CONTRIB to /.bashrc
export FASTJET=$install_dir/fastjet-install
export FASTJET_CONTRIB=$install_dir/fjcontrib-1.046

echo "export FASTJET='$install_dir/fastjet-install'" >> ~/.bashrc
echo "export FASTJET_CONTRIB='$install_dir/fjcontrib-1.046'" >> ~/.bashrc
source ~/.bashrc

# add links to the current directory
cd $current_dir
ln -s $FASTJET/include/fastjet fastjet
