## Velvet Assembler

Velvet is a de novo genomic assembler specially designed for short read sequencing technologies, 
such as Solexa or 454, developed by Daniel Zerbino and Ewan Birney at the European Bioinformatics Institute (EMBL-EBI), near Cambridge, in the United Kingdom.

Velvet currently takes in short read sequences, removes errors then produces high quality unique contigs. It then uses paired-end read and long read information, when available, to retrieve the repeated areas between contigs.

Website: [Velvet Assembler](https://www.ebi.ac.uk/~zerbino/velvet/)

**Mac OSX**:

	brew update
	brew info velvet
	brew install velvet
	
The program is now installed.

**Linux**:

	cd ~/src
	curl -OL https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
	tar zxvf velvet_1.2.10.tgz
	cd velvet_1.2.10
	make

When installed from source as above link the `velveth` and `velvetg` programs to the `bin` folder.

	ln -s ~/src/velvet_1.2.10/velveth ~/bin
	ln -s ~/src/velvet_1.2.10/velvetg ~/bin

Assembling a genome requires running the `velveth` hashing program
followed by the `velvetg` graph assembly.

Test that the programs work with

	velveth
	
and
	
	velvetg
	
	
