MITObim - mitochondrial baiting and iterative mapping
=====================================================



VERSIONS
--------

1.6 (stable - relies on MIRA 3.4.1.1)

1.7 (beta - relies on MIRA 4)


Copyright Christoph Hahn 2012-2014

CONTACT
-------

Christoph Hahn <christoph.hahn@nhm.uio.no>

INTRODUCTION
------------

This document contains instructions on how to use the MITObim pipeline described in Hahn et al. 2013. The full article can be found [here](http://nar.oxfordjournals.org/content/early/2013/05/09/nar.gkt371.full "MITObim full article at NAR"). Kindly cite the article if you are using MITObim in your work. The pipeline is at the moment intended to be used with illumina data, but the use of Iontorrent and 454 data has been enabled in the current versions.


PREREQUISITES
-------------

- GNU utilities
- Perl
- A running version of MIRA 3.4.1.1 (for the use with MITObim 1.6 - download [here](http://sourceforge.net/projects/mira-assembler/files/MIRA/Older%20releases/V3.4.0/)) or MIRA 4 (for the use with MITObim 1.7 - download [here](http://sourceforge.net/projects/mira-assembler/files/MIRA/stable/)) is required. An excellent guide to MIRA is available [here](http://mira-assembler.sourceforge.net/docs/DefinitiveGuideToMIRA.html "The definitive Guide to MIRA").

COMMENT
-------

As I received many requests about enabling MITObim for MIRA 4 I decided to upload a beta version of MITObim 1.7. However, tHe proofreading option is at the moment disabled in this version and will be enabled once I got a chance to thoroughly test its behavior with MIRA 4. If you are planning to use the proofreading functunality please refer to MITObim 1.6 for the time being. SOrry for the inconvenience!


General introduction to MITObim
-------------------------------

The MITObim procedure (mitochondrial baiting and iterative mapping) represents a highly efficient approach to assembling novel mitochondrial genomes of non-model organisms directly from total genomic DNA derived NGS reads. Labor intensive long-range PCR steps prior to sequencing are no longer required. MITObim is capable of reconstructing mitochondrial genomes without the need of a reference genome of the targeted species by relying solely on (a) mitochondrial genome information of more distantly related taxa or (b) short mitochondrial barcoding sequences (seeds), such as the commonly used cytochochrome-oxidase subunit 1 (COI), as a starting reference. 

The script is performing three steps and iteratively repeating them: (i) Deriving reference sequence from previous mapping assembly, (ii) in silico baiting using the newly derived reference (iii) previously fished reads are mapped to the newly derived reference leading to an extension of the reference sequence.
For more details please refer to [Hahn et al. 2013](http://nar.oxfordjournals.org/content/early/2013/05/09/nar.gkt371.full "MITObim full article at NAR"). Detailed examples are demonstrated in the TUTORIALS section below.


TUTORIALS
---------

The following tutorials are designed for users with little Unix and no previous MIRA experience. Tutorials I & II will demonstrate how to recover the complete mitochondrial genome of _Thymallus thymallus_ using the mitochondrial genome of _Salvelinus alpinus_ as a starting reference. Tutorial III achieves the same goal using solely a ~700 bp barcoding sequence as initial seed reference. For a "quick and dirty" exploration of your own data I recommend trying something along the lines of Tutorial II. Tutorial IV (to be added soon) uses a proofreading procedure to specifically reconstruct two mitochondrial genomes from a mixed sample containing genomic reads from two species. FOr convenience I always refer to MITObim.pl in the tutorials - The individual user will have to call the respective version of MITObim (e.g. MITObim_1.6.pl) during trying the tutorials. MITObim 1.7 might finish some of the tutorials with slighly less iterations than described in the tutorials. Dont worry! 


Preparations:

- download/compile MIRA (MIRA 3.4.1.1 for the use with MITObim 1.6 from [here](http://sourceforge.net/projects/mira-assembler/files/MIRA/Older%20releases/V3.4.0/) or MIRA 4 for the use with MITObim 1.7 from [here](http://sourceforge.net/projects/mira-assembler/files/MIRA/stable/)). Precompiled binaries are available for Linux and OSX. Help can be found [here](http://mira-assembler.sourceforge.net/docs/DefinitiveGuideToMIRA.html "Definitive Guide to MIRA"). You ll need to put the directory containing the MIRA executables in your PATH in order to successfully use MITObim.pl. If you can't or won't do that you can also tell MITObim where to find the correct MIRA binaries via the --mirapath option.
- download the MITObim wrapper script and the testdata from Github, e.g. download the entire MITObim repository as zip archive (use the button on the Github page) or use git on the command line (git clone --recursive git://github.com/chrishah/MITObim.git).
- make the MITObim.pl executable (chmod a+x MITObim.pl) and extract the contents of the testdata archives (tar xvfz testdata?.tgz)

Test the wrapper script by doing:

	-bash-4.1$ ~/PATH/TO/MITObim.pl

which should display the usage (NOTE: In MITObim 1.7 the -strain flag has been renamed to -sample):

	usage ./MITObim.pl <parameters>

	parameters:

                -start <int>            iteration to start with, default=1
                -end <int>              iteration to end with, default=1
                -strain <string>        strainname as used in initial MIRA assembly
                -ref <string>           referencename as used in initial MIRA assembly
                -readpool <FILE>        readpool in fastq format
                -maf <FILE>             maf file from previous MIRA assembly

	optional:
                --quick <FILE>          starts process with initial baiting using provided fasta reference
                --kbait <int>           set kmer for baiting stringency (default: 31)
                --denovo                runs MIRA in denovo mode, default: mapping
                --pair                  finds pairs after baiting (relies on /1 and /2 ID convention for read pairs), default: no
                --help                  shows this helpful information
                --clean                 retain only the last 2 iteration directories
                --trim                  trim data (we recommend to trim beforehand and feed MITObim with pre trimmed data)
                --iontor                use data produced by iontorrent (experimental - default is illumina data)
                --454                   use 454 data (experimental - default is illumina data)
                --mirapath <string>     full path to MIRA binaries (only needed if MIRA is not in PATH)
                --proofread             applies proofreading (atm only to be used if starting the process from a single short seed reference)
                --readlength <int>      read length of illumina library, default=150, needed for proofreading
                --insert <int>          insert size of illumina library, default=300, needed for proofreading
	examples:
                ./MITObim.pl -start 1 -end 5 -strain StrainX -ref reference-mt -readpool illumina_readpool.fastq -maf initial_assembly.maf
                ./MITObim.pl -end 10 --quick reference.fasta -strain StrainY -ref reference-mt -readpool illumina_readpool.fastq


The archive testdata1 contains three files:

1. Tthymallus-150bp-300sd50-interleaved.fastq - 6000 simulated illumina reads (read length 150 bp, insert size 300 +- 50 bp) for the mitochondrial genome of _T. thymallus_ as discussed in [Hahn et al](http://nar.oxfordjournals.org/content/early/2013/05/09/nar.gkt371.full "MITObim full article at NAR").

2. Salpinus-mt-genome-NC_000861.fasta - mitochondrial genome of _S. alpinus_ in fasta format downloaded from Genbank (accession NC000861).

3. Tthymallus-COI-partial-HQ961018.fasta - partial COI sequence of _T. thymallus_ (Genbank acc. HQ961018).


TUTORIAL I: reconstruction of a mitochondrial genome using a two step procedure
-------------------------------------------------------------------------------

One strategy for successful MITObim performance is to prepare an initial reference from the conserved regions of a distantly related mitochondrial reference genome using a MIRA mapping assembly in a first step. In a second step the MITObim wrapper script is applied to this newly constructed reference in order to reconstruct the entire mitochondrial genome.

a. Initial mapping assembly using MIRA:

Create a directory to work in and change into this directory:

	-bash-4.1$ mkdir tutorial1
	-bash-4.1$ cd tutorial1

**mapping assembly with MIRA 3.4.1.1:**

MIRA 3.4.1.1 expects certain files to be present in your working directory, which are named according to the name of your MIRA project and certain conventions. Let's say the project name of your choice for the initial mapping assembly is "initial-mapping-testpool-to-Salpinus-mt". One could either simply copy and rename the files accordingly in your working directory or create symbolic links in the current working directory pointing to the actual files (replace "cp" with "ln -s" in the following command). 
  
	-bash-4.1$ cp /PATH/TO/testdata1/Tthymallus-150bp-300sd50-interleaved.fastq initial-mapping-testpool-to-Salpinus-mt_in.solexa.fastq 
	-bash-4.1$ cp /PATH/TO/testdata1/Salpinus-mt-genome-NC_000861.fasta initial-mapping-testpool-to-Salpinus-mt_backbone_in.fasta
  
run MIRA 3.4.1.1:

	-bash-4.1$  mira --project=initial-mapping-testpool-to-Salpinus-mt --job=mapping,genome,accurate,solexa -MI:somrnl=0 -AS:nop=1 -SB:bft=fasta:bbq=30:bsn=Salpinus-mt-genome SOLEXA_SETTINGS -CO:msr=no -SB:dsn=testpool 

or

**mapping assembly with MIRA 4:**

MIRA 4 does not depend on the naming conventions of MIRA 3.4 any more. MIRA 4 uses a so called manifest file in which you can specifiy the paths to your data (see below). Let us for the sake of this example create symbolic links, nevertheless. *NOTE* that we change the extension of the reference fasta file to *.fa. The reason is that MIRA will per default expect quality data for files with the extension *.fasta. The file extension *.fa tells MIRA to continue even without quality data.  

	-bash-4.1$ ln -s /PATH/TO/testdata1/Tthymallus-150bp-300sd50-interleaved.fastq reads.fastq 
	-bash-4.1$ ln -s /PATH/TO/testdata1/Salpinus-mt-genome-NC_000861.fasta reference.fa
	
create a manifest file (named manifest.conf) specifying the parameters for your MIRA assembly (see section 3.4 [here](http://mira-assembler.sourceforge.net/docs/DefinitiveGuideToMIRA.html "Definitive Guide to MIRA") for details) in your favourite text editor or just type (or copy and paste) the following command:

	-bash-4.1$ echo -e "\n#manifest file for basic mapping assembly with illumina data using MIRA 4\n\nproject = initial-mapping-testpool-to-Salpinus-mt\n\njob=genome,mapping,accurate\n\nparameters = -NW:mrnl=0 -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no\n\nreadgroup\nis_reference\ndata = reference.fa\nstrain = Salpinus-mt-genome\n\nreadgroup = reads\ndata = reads.fastq\ntechnology = solexa\nstrain = testpool\n" > manifest.conf

Whatever you do, the manifest.conf file should eventually look as follows:

	-bash-4.1$ head -n 20 manifest.conf

	#manifest file for basic mapping assembly with illumina data using MIRA 4

	project = initial-mapping-testpool-to-Salpinus-mt

	job=genome,mapping,accurate

	parameters = -NW:mrnl=0 -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no

	readgroup
	is_reference
	data = reference.fa
	strain = Salpinus-mt-genome

	readgroup = reads
	data = reads.fastq
	technology = solexa
	strain = testpool

run MIRA 4:

	-bash-4.1$ mira manifest.conf

Look what happened:

	-bash-4.1$ ls -hlrt
	total 2.8M
	-rw-r--r-- 1 chrishah users  17K Oct 21 22:50 initial-mapping-testpool-to-Salpinus-mt_backbone_in.fasta
	-rw-r--r-- 1 chrishah users 2.7M Oct 21 22:51 initial-mapping-testpool-to-Salpinus-mt_in.solexa.fastq
	drwxr-xr-x 6 chrishah users    4 Oct 21 22:51 initial-mapping-testpool-to-Salpinus-mt_assembly

The successfull MIRA instance created a directory "initial-mapping-testpool-to-Salpinus-mt_assembly", containing four directories:

	-bash-4.1$ ls -hlrt initial-mapping-testpool-to-Salpinus-mt_assembly/
	total 2.0K
	drwxr-xr-x 2 chrishah users  2 Oct 21 22:51 initial-mapping-testpool-to-Salpinus-mt_d_chkpt
	drwxr-xr-x 2 chrishah users 15 Oct 21 22:51 initial-mapping-testpool-to-Salpinus-mt_d_tmp
	drwxr-xr-x 2 chrishah users  9 Oct 21 22:51 initial-mapping-testpool-to-Salpinus-mt_d_results
	drwxr-xr-x 2 chrishah users  8 Oct 21 22:51 initial-mapping-testpool-to-Salpinus-mt_d_info
 
The newly constructed reference is contained in the file "initial-mapping-testpool-to-Salpinus-mt_out.maf" in the "initial-mapping-testpool-to-Salpinus-mt_d_results" directory.

b. Baiting and iterative mapping using the MITObim.pl script: 

Run the wrapper script as follows (standard output will be written into the file log - *approximate runtime: 2 min*):

for **MITObim 1.6**:

	-bash-4.1$ /PATH/TO/MITObim.pl -start 1 -end 10 -strain testpool -ref Salpinus_mt_genome -readpool initial-mapping-testpool-to-Salpinus-mt_in.solexa.fastq -maf initial-mapping-testpool-to-Salpinus-mt_assembly/initial-mapping-testpool-to-Salpinus-mt_d_results/initial-mapping-testpool-to-Salpinus-mt_out.maf &> log

for **MITObim 1.7**:

	-bash-4.1$ /PATH/TO/MITObim.pl -start 1 -end 10 -sample testpool -ref Salpinus_mt_genome -readpool reads.fastq -maf initial-mapping-testpool-to-Salpinus-mt_assembly/initial-mapping-testpool-to-Salpinus-mt_d_results/initial-mapping-testpool-to-Salpinus-mt_out.maf &> log

NOTE: 
- The strain/sample name needs to be the same as used in the initial MIRA assembly.

After the process has finished looking into the log file will yield something like:

	-bash-4.1$ tail log

	End of assembly process, thank you for using MIRA.

	readpool contains 6000 reads
	contig length: 16667

	MITObim has reached a stationary read number after 8 iterations!!

Allthough you have asked for 10 iterations MITObim stopped the process after fewer iterations because it has detected a stationary readnumber. It will have created a number of iteration* directories. Each of which is organized like the working directory above (initial mapping assembly), containing a MIRA assembly directory plus the referernce and readpool used in the respective iteration:

	-bash-4.1$ ls -hlrt iteration1/
	total 2.3M
	-rw-r--r-- 1 chrishah users  17K Oct 21 23:45 testpool-Salpinus_mt_genome_backbone_in.fasta
	-rw-r--r-- 1 chrishah users 383K Oct 21 23:45 hashstat.bin
	-rw-r--r-- 1 chrishah users 1.9M Oct 21 23:45 testpool-Salpinus_mt_genome_in.solexa.fastq
	drwxr-xr-x 6 chrishah users    4 Oct 21 23:45 testpool-Salpinus_mt_genome_assembly

A fasta file containing the complete mitochondondrial genome of _T. thymallus_ can be found in the final iteration directory. See its contents by doing e.g.:

	-bash-4.1$ less iteration8/testpool-Salpinus_mt_genome_assembly/testpool-Salpinus_mt_genome_d_results/testpool-Salpinus_mt_genome_out.unpadded.fasta

**Congratulations!!!**

*IMPORTANT NOTE*: For real data this strategy might require substantial computational resources depending on the size of the illumina readpool (Memory requirements can be estimated using the miramem program shipping with MIRA). This increased memory consumption can be bypassed by limitating the readpool e.g. via an initial in-silico baiting step using mirabait (in-silico baiting program shipping with MIRA). While this strategy will impair the initial sensitivity of the approach and might therefore require more iterations to finish, it works reasonably well if a not too distantly related reference is available. This strategy can be performed by using the --quick flag, together with providing a reference sequence in fasta format. A detailed example can be found in TUTORIAL II below.

TUTORIAL II - direct reconstruction without prior mapping assembly using the --quick option
-------------------------------------------------------------------------------------------
This TUTORIAL illustrates the "quick & dirty" strategy that I usually use in a first test. To finish a mitochondrial genome it usually takes more iterations than TUTORIAL I above because the initial mapping assembly is less thorough. The **big** advantage is that the whole process can be run on a standard Desktop computer due to the substantial reduction in the number of reads to be dealt with already in the first iteration.
Run the MITObim.pl script with the --quick option, providing a reference in fasta format (*approximate runtime: 4 min*):

	-bash-4.1$ mkdir tutorial2
	-bash-4.1$ cd tutorial2

for **MITObim 1.6**:

	-bash-4.1$ /PATH/TO/MITObim.pl -start 1 -end 30 -strain testpool -ref Salpinus_mt_genome -readpool ~/PATH/TO/testdata1/Tthymallus-150bp-300sd50-interleaved.fastq --quick ~/PATH/TO/testdata1/Salpinus-mt-genome-NC_000861.fasta &> log

for **MITObim 1.7**:

	-bash-4.1$ /PATH/TO/MITObim.pl -start 1 -end 30 -sample testpool -ref Salpinus_mt_genome -readpool ~/PATH/TO/testdata1/Tthymallus-150bp-300sd50-interleaved.fastq --quick ~/PATH/TO/testdata1/Salpinus-mt-genome-NC_000861.fasta &> log

You will find that with this approach MITObim will reconstruct the mitochondrial genome and reach a stationary number of mitochondrial reads only after 15 (or less) iterations. The result should nevertheless be equal to that obtained in TUTORIAL I and we were able to bypass the increased memory requirements of the initial MIRA mapping assembly.

**CONGRATULATIONS!!** ..but, there is more..

TUTORIAL III - reconstructing mt genomes from mt barcode seeds
--------------------------------------------------------------

This tutorial reconstructs the mt genome of _T. thymallus_ solely using a partial mitochondrial COI sequence as starting seed. *NOTE* that we also use the new --clean option which tells MITObim to always only keep the latest two iteration directories to save space (*approximate runtime: 20 min*):

	-bash-1.4$ mkdir tutorial3
	-bash-1.4$ cd tutorial3

for **MITObim 1.6**:

	-bash-1.4$ ~/PATH/TO/MITObim.pl -strain testpool -ref Tthymallus-COI -readpool ~/PATH/TO/testdata1/Tthymallus-150bp-300sd50-interleaved.fastq --quick ~/PATH/TO/testdata1/Tthymallus-COI-partial-HQ961018.fasta -end 100 --noshow --clean &> log

for **MITObim 1.7**:

	-bash-1.4$ ~/PATH/TO/MITObim.pl -sample testpool -ref Tthymallus-COI -readpool ~/PATH/TO/testdata1/Tthymallus-150bp-300sd50-interleaved.fastq --quick ~/PATH/TO/testdata1/Tthymallus-COI-partial-HQ961018.fasta -end 100 --clean &> log

MITObim reconstructs the mitchondrial genome in 82 (or less) iterations.

For "well behaved" datasets the standard mapping assembly can be substituted by a _de novo_ assembly (--denovo flag). Utilizing read pair information (--paired flag) can further speed up the reconstruction if run in _de novo_ mode. This however will not decrease the nuber of necessary iterations in standard mapping mode.
Test this strategy, like so (*approximate runtime: 10 min*):

	-bash-4.1$ mkdir tutorial3-denovo
	-bash-4.1$ cd tutorial3-denovo

for **MITObim 1.6**:

	-bash-4.1$ ~/PATH/TO/MITObim.pl -strain testpool -ref Tthymallus-COI -readpool ~/PATH/TO/testdata1/Tthymallus-150bp-300sd50-interleaved.fastq --quick ~/PATH/TO/testdata1/Tthymallus-COI-partial-HQ961018.fasta -end 50 --noshow --denovo --paired --clean &> log

for **MITObim 1.7**:

	-bash-4.1$ ~/PATH/TO/MITObim.pl -sample testpool -ref Tthymallus-COI -readpool ~/PATH/TO/testdata1/Tthymallus-150bp-300sd50-interleaved.fastq --quick ~/PATH/TO/testdata1/Tthymallus-COI-partial-HQ961018.fasta -end 50 --denovo --paired --clean &> log

This strategy reconstructs the correct mitochondrial genome in only 31 (or less) iterations.


