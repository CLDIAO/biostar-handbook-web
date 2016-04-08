## Bedops

BEDOPS is an open-source command-line toolkit that performs highly efficient
and scalable Boolean and other set operations, statistical calculations, archiving,
conversion and other management of genomic data of arbitrary scale. Tasks can be
easily split by chromosome for distributing whole-genome analyses across a
computational cluster.

Website: http://bedops.readthedocs.org/en/latest/

**Mac OS X**:

    brew update
    brew info bedops
    brew install bedops

**Linux**:

    wget https://github.com/bedops/bedops/archive/master.zip && unzip master.zip
    cd bedops-master
    make && make install
    cp bin/* /usr/local/bin

**Linux and Mac**

See the BEDOPS release page: https://github.com/bedops/bedops/releases

Test installation by running:

	bedops
