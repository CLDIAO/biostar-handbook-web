## BLAST

BLAST stands for **Basic Local Alignment Search Tool**
and is an algorithm for searching a sequence (called a query)
against very large sequence databases (called a target).

NCBI has is a [web interface][web-blast] for BLAST and provides
[downloadable blast][local-blast] programs.

Installation on **Mac OSX**:

	brew update
	brew info blast
	brew install blast

Installation on **Ubuntu and Windows Bash**

    sudo apt-get install -y ncbi-blast+

Test the version of blast that you have:

    blastn -version

It will print:

    blastn: 2.2.28+
    Package: blast 2.2.28, build Jun  3 2013 11:17:14

Often there are newer versions for the BLAST software that have added
functionality. If we need these new features we must install blast manually.

Investigate the version of the [latest release][local-blast]
and adapt the names of the programs accordingly. These
instructions below will install version `BLAST+ 2.3.0`

The URLs below will change, make sure to set them correctly

**Mac OSX**:

    URL=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.3.0+-universal-macosx.tar.gz

**Ubuntu and Windows Bash**:

    URL=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz

Then install in the following way:

    cd ~/src
    curl -O  $URL
    tar zxvf ncbi-blast-*.tar.gz
    export PATH=~/src/ncbi-blast-2.3.0+/bin:$PATH


**Windows**:

Blast can be easily installed on Windows system.

* Visit [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* Download the file `ncbi-blast-2.2.31+-win64.exe`
* Double-click and install

**Note** that installing into the default path can make
blast fail in some circumstances because this path includes
a directory named `Program Files` and the space in the name of the
directory can cause various problems later.
It is best if the installation is performed into the cygwin directory:

    C:\cygwin64\home\your-user-name\src\

New Cygwin terminals will be able to run all blast tools.

[web-blast]: http://blast.ncbi.nlm.nih.gov/Blast.cgi
[local-blast]: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
