## Subread and featureCounts

Aligner and feature counter

Website: http://bioinf.wehi.edu.au/subread-package/

Obtain the source code:

	cd ~/src
	curl -OL http://sourceforge.net/projects/subread/files/subread-1.5.0/subread-1.5.0-source.tar.gz
	tar zxv subread-1.5.0-source.tar.gz
	cd subread-1.5.0-source/src

**Mac OSX**:

	make -f Makefile.MacOS

**Linux** :

    make -f Makefile.Linux

Link the executables into the `~bin` folder:

	ln -fs ~/src/subread-1.5.0-source/bin/featureCounts ~/bin
	ln -fs ~/src/subread-1.5.0-source/bin/subread-align ~/bin
	ln -fs ~/src/subread-1.5.0-source/bin/subread-buildindex ~/bin

Check that the programs work

    featureCounts
    subread-align

