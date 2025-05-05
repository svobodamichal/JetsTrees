## FASTJET install 


This framework requires the installation of FastJet + FastJet contrib.
FastJet can be found at: http://fastjet.fr/all-releases.html
* download FastJet 3.0 +

FastJet contrib can be found at:
```
http://fastjet.hepforge.org/contrib/
```
RCAS & PDSF both require 32-bit build to fit into other software framework.  To force this install on the 64-bit systems, please follow the below prescription to avoid headaches.

FASTJET installation example:
```bash
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
```

FastJet contrib installation example:
```bash
    cd  $install_dir
    wget http://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.046.tar.gz
    tar zxvf fjcontrib-1.046.tar.gz
    cd fjcontrib-1.046
    ./configure --fastjet-config=$PWD/../fastjet-3.4.2/fastjet-config --prefix=$PWD/../fastjet-install CXXFLAGS="-m32 -fPIC -fno-inline" CFLAGS="-m32 -fPIC -fno-inline" LDFLAGS="-m32"
    make
    make install
    make fragile-shared
    make fragile-shared-install
```
unset flags after:
```bash
  unset CXXFLAGS CFLAGS LDFLAGS
```
You will additionally need to set FastJet environmental variables in your startup file (mine is `.bashrc`, and the below is an example):
```bash
export FASTJET='/gpfs01/star/pwg/$USER/install/fastjet-install'
export FASTJET_CONTRIB='/gpfs01/star/pwg/$USER/install/fjcontrib-1.046'
echo "export FASTJET='$install_dir/fastjet-install'" >> ~/.bashrc
echo "export FASTJET_CONTRIB='$install_dir/fjcontrib-1.046'" >> ~/.bashrc
source ~/.bashrc
```
and soft link the include paths to your working directory which contains your `StRoot `directory and where you run your `macro.C` file from:

```bash
cd $current_dir
ln -s $FASTJET/include/siscone siscone
ln -s $FASTJET/include/fastjet fastjet
```

# Embedding