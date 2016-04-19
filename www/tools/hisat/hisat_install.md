## HISAT2 aligner

HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads 
(both DNA and RNA) against the general human population (as well as against a single reference genome). 

Website: [HiSat Website](https://ccb.jhu.edu/software/hisat2/index.shtml)

**Mac OSX**:

	brew update
	brew info hisat2
	brew install hisat2
	
The program is now installed.

**Linux**

Find and download the Linux version (also applies as an alternative on Mac OSX)

	cd ~/src
	curl -OL ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.3-beta-Linux_x86_64.zip
	unzip hisat2-2.0.3-beta-Linux_x86_64.zip
	
The **Linux** link the executables:

	ln -s ~/src/hisat2-2.0.3-beta/hisat2-build ~/bin
	ln -s ~/src/hisat2-2.0.3-beta/hisat2 ~/bin
	
Test that the programs work with

	hisat2
	
	
	
