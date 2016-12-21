#! /usr/bin/perl
#
# MITObim - mitochondrial baiting and iterative mapping
# wrapper script version 1.8
# Author: Christoph Hahn, 2012-2016
# christoph.hahn@nhm.uio.no

use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Copy;
use List::Util qw< min max >;
use POSIX qw(strftime);
use POSIX qw(ceil);
use File::Path 'rmtree';

my $startiteration = 1;
my $enditeration = 1;

my ($quick, $noshow, $help, $strainname, $paired, $mode, $refname, $readpool, $maf, $proofreading, $readlength, $insertsize, $MM, $trim, $trimoverhang, $k_bait, $clean, $clean_interval, $min_contig_cov) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 31, 0, 2, 0);
my ($miramode, $key, $val, $exit, $current_contiglength, $current_number_of_reads, $current_number_of_contigs);
my $splitting = 0;
my $platform = "solexa";
my $platform_settings;# = "SOLEXA";
my $shme = "";
my $trim_off = "";
my $redirect_temp = "";
my $NFS_warn_only = "";
my ($mirapath, $mira, $miraconvert, $mirabait) = ("", "mira", "miraconvert", "mirabait");
my (@reads, @output, @path, @current_contig_stats, @contiglengths, @number_of_reads, @number_of_contigs);
my %hash;
my $PROGRAM = "\nMITObim - mitochondrial baiting and iterative mapping\n";
my $VERSION = "version 1.8\n";
my $AUTHOR = "author: Christoph Hahn, (c) 2012-2016\n";
my $cite = "\nif you found MITObim useful, please cite:
Hahn C, Bachmann L and Chevreux B. (2013) Reconstructing mitochondrial genomes directly from genomic next-generation sequencing reads -
a baiting and iterative mapping approach. Nucl. Acids Res. 41(13):e129. doi: 10.1093/nar/gkt371\n\n";
my $USAGE = 	"\nusage: ./MITObim.pl <parameters>
	 	\nparameters:
		-start <int>		iteration to start with, default=1
		-end <int>		iteration to end with, default=1
		-sample <string>	sampleID as used in initial MIRA assembly
		-ref <string>		referenceID as used in initial MIRA assembly
		-readpool <FILE>	readpool in fastq format (*.gz is also allowed)
		-maf <FILE>		maf file from previous MIRA assembly
		\noptional:
		--quick <FILE>		starts process with initial baiting using provided fasta reference
		--kbait <int>		set kmer for baiting stringency (default: 31)
		--platform		specify sequencing platform (default: 'solexa'; other options: 'iontor', '454', 'pacbio')
		--denovo		runs MIRA in denovo mode (default: mapping)
		--pair			finds pairs after baiting (relies on /1 and /2 header convention for read pairs) (default: no)
		--verbose		show detailed output of MIRA modules (default: no)
		--split			split reference at positions with more than 5N (default: no)
		--help			shows this helpful information
		--clean                 retain only the last 2 iteration directories (default: no)
		--trimreads		trim data (default: no; we recommend to trim beforehand and feed MITObim with pre trimmed data)
		--trimoverhang		trim overhang up- and downstream of reference (default: no)
		--missmatch <int>	number of allowed missmatches in mapping - only for illumina data (default: 15% of avg. read length)
		--min_cov <int>		minimum average coverage of contigs to be retained (default: off)
		--mirapath <string>     full path to MIRA binaries (only needed if MIRA is not in PATH)
		--redirect_tmp		redirect temporary output to this location (useful in case you are running MITObim on an NFS mount)
		--NFS_warn_only		allow MIRA to run on NFS mount without aborting -  warn only (expert option - see MIRA documentation 'check_nfs')
		\nexamples:
		./MITObim.pl -start 1 -end 5 -sample StrainX -ref reference-mt -readpool illumina_readpool.fastq -maf initial_assembly.maf
		./MITObim.pl -end 10 --quick reference.fasta -sample StrainY -ref reference-mt -readpool illumina_readpool.fastq\n\n";
#		--proofread		applies proofreading (atm only to be used if starting the process from a single short seed reference)
#		--readlength <int>	read length of illumina library, default=150, relevant only for proofreading
#		--insert <int>		insert size of illumina library, default=300, relevant only for proofreading

my $command = $0;
for (@ARGV){
        $command .= " $_";
}

GetOptions (	"start=i" => \$startiteration,
		"end=i" => \$enditeration,
		"kbait=i" => \$k_bait,
		"quick=s" => \$quick,
		"verbose!" => \$noshow,
		"sample=s" => \$strainname,
		"paired" => \$paired,
		"denovo" => \$mode,
		"ref=s" => \$refname,
		"readpool=s" => \$readpool,
		"help!" => \$help,
		"clean!" => \$clean,
		"mirapath=s" => \$mirapath,
		"maf=s" => \$maf,
#		"proofreading!" => \$proofreading,
		"trimreads!" => \$trim,
		"trimoverhang!" => \$trimoverhang,
		"missmatch=i" => \$MM,
		"platform=s" => \$platform,
#		"readlength=i" => \$readlength,
#		"insertsize=i" => \$insertsize,
		"split!"	=>	\$splitting,
		"min_cov=i"	=>	\$min_contig_cov,
		"redirect_tmp=s" =>	\$redirect_temp,
		"NFS_warn_only!" => \$NFS_warn_only) or die "Incorrect usage!\n$USAGE";


print $PROGRAM; 
print $VERSION; 
print $AUTHOR; 

print $USAGE and exit if $help;
print $USAGE and exit if !$readpool;
unless ($quick){
        print $USAGE and exit if !$maf;
}
print $USAGE and exit if !$refname;

$readpool=abs_path($readpool);
unless (-e $readpool){
	print "Cant find the readpool. Is the path correct?\n";
}
if ($maf){
	$maf=abs_path($maf);
	unless (-e $maf){
		print "Cant find *.maf file. Is the path correct?\n";
	}
}

($platform, $platform_settings) = &set_platform($platform);

if ($quick){
	$quick=abs_path($quick);
        unless (-e $quick){
		print "quick option selected but is the path to the file correct?\n";
		exit 1;
	}
	print "\nquick option selected! -maf option will be ignored (if given)\n";
	$maf = 0;
	$startiteration = 0;
}
print $USAGE and exit if ($startiteration > $enditeration);
unless (((-e $maf)||($quick)) && (-e $readpool)){
        print "\nAre readpool AND maf/reference files there?\n";
        exit 1;
}
if ($mirapath){
	if (-e "$mirapath/mira"){
		print "found executables in the path specified by the user - good!\n";
		$mira = "$mirapath/mira";
		$miraconvert = "$mirapath/miraconvert";
		$mirabait = "$mirapath/mirabait";
	}else{
		print "somethings wrong with the path to mira.\n";
		exit 1;
	}
}


##if not given otherwise, readlength and insertsize are set to default. automatic readlength and insertsize detection will be implemented in time.
#if (!$readlength){
#	$readlength = 150;
#}
#if (!$insertsize){
#	$insertsize = 300;
#}
if (!$mode){
	$miramode = "mapping";
}else {
	$miramode = "denovo";
}

if (!$trim){
	$trim_off = "\"--noclipping -CL:pec=no\"";
}


print "\nFull command run:\n$command\n";
print "\nAll paramters seem to make sense:\n";
print "startiteration: $startiteration\n";
print "enditeration: $enditeration\n";
print "sample: $strainname\n";
print "refname: $refname\n";
print "readpool: $readpool\n";
print "maf: $maf\n";
print "quick: $quick\n";
print "paired: $paired (off=0, on=1)\n";
print "assembly mode: $mode (mapping=0, denovo=1)\n";
print "verbose: $noshow (off=0, on=1)\n";
print "split: $splitting (off=0, on=1)\n";
print "minimum avg. coverage: $min_contig_cov (off=0)\n";
print "clean: $clean (off=0, on=1)\n";
print "trim reads: $trim (off=0, on=1)\n";
print "trim overhang: $trimoverhang (no=0, yes=1)\n";
print "platform: $platform\n";
print "kmer baiting: $k_bait\n";
#print "proofread: $proofreading (off=0, on=1)\n";

if ($proofreading){
	print "\nproofreading is not yet enabled in this beta version of MITObim 1.8 - it is currently being optimized for MIRA4. \nplease refer to MITObim 1.6 if you wish to use this option. Sorry for the inconvenience!\n\n";
	exit;

	print "proofreading: on\n";
	print "readlength: $readlength\n";
	print "insertsize: $insertsize\n";
	$MM = 0;
	print "number of allowed missmatches in proofreading assembly: $MM\n";
	$shme = "-AL:shme=$MM";
}elsif ((!$proofreading) && (!$mode) && ($platform eq "solexa")){
	if ($MM == -1){
		print "number of missmatches in mapping assembly: default (15% of average read length loaded)\n";
		$shme = "";
	}else {
		print "number of missmatches in mapping assembly: $MM\n";
		$shme = "-AL:shme=$MM";
	}
	print "proofreading: off\n";
}

if (!$trimoverhang){
	$trimoverhang = "-SB:tor=no";
}else {
	$trimoverhang = "-SB:tor=yes";
}

if ($redirect_temp) {
	$redirect_temp="-DI:trt=$redirect_temp"
}

print "\nStarting MITObim \n";

my @iteration = ($startiteration .. $enditeration);
foreach (@iteration){
	chomp;
	my $currentiteration = $_;
	mkdir "iteration$currentiteration" or die "MITObim will not overwrite an existing directory: iteration$currentiteration\n";
	chdir "iteration$currentiteration" or die $!;
	print "\n==============\n";
	print " ITERATION $currentiteration\n";
	print "==============\n";
	print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
#	if (($proofreading) && ($currentiteration != 0)){
	if ($proofreading){
		$shme = "-AL:shme=0";
	} 
		
	if ($maf){
		print "\nrecover backbone by running miraconvert on maf file\n\n";
		if ($currentiteration<2){
			@output= qx($miraconvert -f maf -t fasta -A "$platform_settings -CO:fnicpst=yes" $maf tmp);
		}else{
			@output= qx($miraconvert -f maf -t fasta -y $min_contig_cov -A "$platform_settings -CO:fnicpst=yes" $maf tmp);
		}
		$exit = $? >> 8;
		unless (!$noshow){
			print "@output\n";
		}
		unless ($exit == 0){
			if (!$noshow){
				print "@output\n";
			}
			print "\nmiraconvert seems to have failed - see detailed output above\n";
			exit;
		}
#		if ( ((($mode) && ($currentiteration > 1)) && (!$quick)) || ((($mode) && ($currentiteration >= 1)) && ($quick)) ){
#			open(FH1,"<tmp_$strainname.unpadded.fasta") or die "$!";
#		}else{
		open(FH1,"<tmp_$strainname.unpadded.fasta") or die "$!\nIs the sampleID identical as in the initial MIRA assembly?\n";
#		}
		open(FH2,">backbone_it$currentiteration\_initial_$refname.fna") or die $!;
		while (<FH1>) {
			$_ =~ s/x/N/g;
			print FH2 $_; 
		}
		close(FH1);
		close(FH2);
		unlink glob ("tmp*");
	}
	MIRABAIT:
	unless ($maf){
		print "\nquick option baits reads from provided reference in iteration 0\n";
		copy("$quick", "backbone_it$currentiteration\_initial_$refname.fna") or die "copy failed: $!";		
	}
	&check_ref_length("backbone_it$currentiteration\_initial_$refname.fna","temp_baitfile.fasta",29800,$k_bait);

	if ($splitting){
		&remove_unmapped_contigs("temp_baitfile_backup.fasta","temp_baitfile.fasta");
		if (-e "temp_baitfile_backup.fasta"){	#this file has only been created if splitting contigs was necessary
			copy "temp_baitfile.fasta","backbone_it$currentiteration\_initial_$refname.fna";
		}
	}
	print "\nfishing readpool using mirabait (k = $k_bait)\n\n";
	
#	@output = (`mirabait -k $k_bait -n 1 temp_baitfile.fasta $readpool $strainname-$refname\_in.$platform 2>&1`);
	@output = qx($mirabait -k $k_bait -n 1 temp_baitfile.fasta $readpool $strainname-readpool-it$currentiteration);
	$exit = $? >> 8;
	unless (!$noshow){
		print "@output\n";
	}
	if (!-s "$strainname-readpool-it$currentiteration.fastq"){
		print "\nyour readpool does not contain any reads with reasonable match (k = $k_bait) to your reference - Maybe you ll want to different settings or even a different reference?\n\n";
		exit;
	}
	unless ($exit == 0){
		if (!$noshow){
			print "@output\n";
		}
	        print "\nmirabait seems to have failed - see detailed output above\n";
	        exit;
	}
	
	FINDPAIRS:
	
	unless (!$paired){
		print "\nfinding pairs for baited reads\n\n";
		copy ("$strainname-readpool-it$currentiteration.fastq", "$strainname-readpool-it$currentiteration-se.fastq") or die "copy failed: $!";
		open(FH1,"<$strainname-readpool-it$currentiteration.fastq") or die $!;
		open(FH2,">list");
		my $index=1;
		while (<FH1>) {
			if ($index % 8 ==1 || $index % 8 ==5) {
				chomp;
				$_ =~ s/@//g;
				($key, $val) = split /\//;
				$hash{$key} .= exists $hash{$key} ? ",$val" : $val;
			}
			$index++;
		}

		for (keys %hash){
			$_ =~ s/$/\/1/g;
			print FH2 "$_\n";
			$_ =~ s/1$/2/g;
	        	print FH2 "$_\n";
		}
		close(FH1);
		close(FH2);
#		@output = (`miraconvert -f fastq -t fastq -n list $readpool $strainname-$refname\_in.$platform 2>&1`);
		@output = qx($miraconvert -f fastq -t fastq -n list $readpool $strainname-readpool-it$currentiteration );
		$exit = $? >> 8;
		unless (!$noshow){
			print "@output\n";
		}
		unless ($exit == 0){
			if (!$noshow){
				print "@output\n";
			}
	        	print "\nmiraconvert seems to have failed - see detailed output above\n";
	        	exit;
		}
	}
	unlink("list");
	
	MIRA:
	print "\nrunning $miramode assembly using MIRA\n\n";
	&create_manifest($currentiteration,$strainname,$refname,$miramode,$trim_off,$platform_settings,$shme,$paired,$trimoverhang,"$strainname-readpool-it$currentiteration.fastq","backbone_it$currentiteration\_initial_$refname.fna", $redirect_temp, $NFS_warn_only);
	@output = qx($mira manifest.conf ); 

	$exit = $? >> 8;
	unless (!$noshow){
		print "@output\n";
	}
	unless ($exit == 0){
		if (!$noshow){
			print "@output\n";
		}
		print "\nMIRA seems to have failed - see detailed output above\n";
		exit;
	}
	
	if ($redirect_temp) {
		my $tmp_link = readlink("$strainname-$refname\_assembly/$strainname-$refname\_d_tmp");
		rmtree ("$strainname-$refname\_assembly/$strainname-$refname\_d_tmp") or die $!;
		mkdir "$strainname-$refname\_assembly/$strainname-$refname\_d_tmp/" or die $1; 
		my @files = glob("$tmp_link/*");			
		for (@files) {
			move ("$_", "$strainname-$refname\_assembly/$strainname-$refname\_d_tmp/") or die $!;
		}
		rmtree ("$tmp_link") or die $!;
	}
	
	@path = abs_path;
        push (@path, "/$strainname-$refname\_assembly/$strainname-$refname\_d_results/$strainname-$refname\_out.maf");
        $maf = join("",@path);
        unless (-e $maf){
                print "maf file is not there \n";
                exit;
        }
	
#	$current_contiglength = &get_contig_length("$strainname-$refname\_assembly/$strainname-$refname\_d_info/$strainname-$refname\_info_contigstats.txt");
#	$current_number_of_reads = (&get_number_of_reads("$strainname-$refname\_assembly/$strainname-$refname\_d_info/$strainname-$refname\_info_contigstats.txt") - 1);

	@current_contig_stats = &get_contig_stats("$strainname-$refname\_assembly/$strainname-$refname\_d_info/$strainname-$refname\_info_contigstats.txt");
	if (((scalar @current_contig_stats > 3) || ($current_contig_stats[0] > 1)) && ($proofreading)) {
		print "assembly consists of more than one contigs - this is atm not permitted in proofreading mode. Sorry!\n\n";
		exit 1;
	}

	PROOFREAD:
#	if (($proofreading) && ($currentiteration >= 1)){
	if ($proofreading){
		print "proofreading option is currently disabled in this beta version of MITObim 1.7 - sorry for the inconvenience!\n\n";
		exit;
	}

	$current_number_of_contigs = shift @current_contig_stats;
	$current_number_of_reads = shift @current_contig_stats;
	if (!$mode){ #in mapping assemblies the reference is counted as one read
		$current_number_of_reads -= $current_number_of_contigs;
	}
	push (@number_of_contigs, $current_number_of_contigs);
	push (@number_of_reads, $current_number_of_reads);
	print "readpool contains $current_number_of_reads reads\n";
	print "assembly contains $current_number_of_contigs contig(s)\n";
	if (scalar @current_contig_stats == 1){
		print "contig length: $current_contig_stats[0]\n";
	}elsif (scalar @current_contig_stats == 3){
		print "min contig length: ".$current_contig_stats[0]." bp\nmax contig length: ".$current_contig_stats[1]." bp\navg contig length: ".sprintf("%.0f", $current_contig_stats[2])." bp\n";
		print "find details on individual contigs in: ". abs_path . "/$strainname-$refname\_assembly/$strainname-$refname\_d_info/$strainname-$refname\_info_contigstats.txt\n";
#		printf("%1f",$current_contig_stats[2])."\n";
	}else {
		print "somethings wrong with your contig stats. Sorry!\n";
		exit 1;
	}

        if ($clean){
                &clean($clean_interval, $currentiteration);
        }

	if ($number_of_reads[-2]){
		if ($number_of_reads[-2] >= $number_of_reads[-1]){
			print "\nMITObim has reached a stationary read number after $currentiteration iterations!!\n$cite";
			print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
			exit;
		}
	}
	chdir ".." or die "Failed to go to parent directory: $!";
}
print "\nsuccessfully completed $enditeration iterations with MITObim! " . strftime("%b %e %H:%M:%S", localtime) . "\n$cite";

#
#
###SUBROUTINES
#
#
#

sub set_platform{
	my @platforms = ('solexa', '454', 'iontorrent', 'pacbio');
	my %pf_settings = ("solexa" => "SOLEXA_SETTINGS", "454" => "454_SETTINGS", "iontorrent" => "IONTOR_SETTINGS", "pacbio" => "PCBIOHQ_SETTINGS -CO:mrpg=5");
	my $user_setting = shift;
	my @user_choice;
	for (@platforms){
		if ( $_ =~ /^$user_setting/ ){
			push(@user_choice, $_);
		}
	}

	if (scalar(@user_choice) > 1){
		print "\nYour platform choice was ambiguous: $user_setting could be: @user_choice - Please try again\n";
		exit;

	} elsif (scalar(@user_choice) == 0){
		print "\nPlease specify a valid sequencing platform - not specifying anything will default into: $platforms[0]\n\navailable options are:\n@platforms\n";
		exit;
	}
	
	return ($user_choice[0], $pf_settings{$user_choice[0]})
}

sub get_contig_stats{
        my $contig = $_[0];
        my @array;
        my @contiglength;
        my (@readnumber, @stats);
	my $readssum = 0;
        open (CONTIGSTATS,"<$contig") or die $!;
        while (<CONTIGSTATS>){
                unless ($_ =~ /#/){
                        @array = split /\t/;
                        push (@contiglength, $array[1]);
                        push (@readnumber, $array[3]);
                }
        }
        close (CONTIGSTATS);
	if (scalar @readnumber == 1){
		push (@stats, (scalar @readnumber, $readnumber[0], $contiglength[0])); #@stats contains: number of contigs, total number of reads used to build the contigs, length of contig
	}
	elsif (scalar @readnumber > 1){
		$readssum += $_ for @readnumber;
		my $minlength = min @contiglength;
		my $maxlength =	max @contiglength;		
		my @avglength = &standard_deviation(@contiglength);
		push (@stats, (scalar @readnumber, $readssum, $minlength, $maxlength, $avglength[0])); #@stats contains: number of contigs, total number of reads used to build the contigs, minimal, maximal, avg length of contigs
	}
        return @stats;
}


#sub get_contig_length{
#	my $contig = $_[0];
#	my @contiglength;
#	open (CONTIGSTATS,"<$contig") or die $!;
#	while (<CONTIGSTATS>){
#		unless ($_ =~ /#/){
#			@contigslength = split /\t/;
#		} 
#	}
#	close (CONTIGSTATS);		
#	return $contiglength[1];
#}

#sub get_number_of_reads{
#	my $contig = $_[0];
#	my @contigstats;
#	open (CONTIGSTATS,"<$contig") or die $!;
#	while (<CONTIGSTATS>){
#		unless ($_ =~ /#/){
#			@contigstats = split /\t/;
#		} 
#	}
#	close (CONTIGSTATS);		
#	return $contigstats[3];
#}

sub proofread {
        my $zero_MM = $_[0];
        my $readtaglist_FH = $_[1];
	my $contiglength = $_[2];
	my $elevated_cov_lower_start = $_[3];
	my $elevated_cov_lower_end = $_[4];
	my $elevated_cov_upper_start = $_[5];
	my $elevated_cov_upper_end = $_[6];
        my $lower_limit = $_[7];
#        my $lower_limit = 200;
	my $lower_main_limit = $_[8];
#	my $lower_main_limit = 500;
        my $upper_limit = $contiglength - $lower_limit;
        my $upper_main_limit = $contiglength - $lower_main_limit;
	my @readtaglist;
        my $ref;
        my $junk;
        my $current_id;
        my %count;
        my @readid_good;
        my @readid_bad;
        my @reads;
        my @readlist;
        my @readlist_good;
        my @readlist_proofread;
        my @total_readlist;
        my @singleton;
        my @singletons;
        my @taglist;
        my @taglist_line;
        my @readtaglist_lower;
        my @readtaglist_upper;
        my @read_ids_lower;
        my @read_ids_all;
        my @read_ids_upper;
        my %ids =();
        my @unsorted;
        my $min;
        my $max;
        my $tag;


	print "\nlower limit: $lower_limit\n";
        print "upper limit: $upper_limit\n";
	print "lower main limit: $lower_main_limit\n";
	print "upper main limit: $upper_main_limit\n\n";
        open (TAGLIST,"<$readtaglist_FH") or die $!;
        while (<TAGLIST>){
                push (@readtaglist, "$_");
        }
        close (TAGLIST);

        for (@readtaglist){
                @taglist_line = split /\t/;
                unless ($taglist_line[0] =~ /#/){
                        $ref = join ("\t", $taglist_line[0], $taglist_line[6]);
                        push (@read_ids_all, $ref);

			if (($taglist_line[1] <= $lower_limit) || ($taglist_line[2] <= $lower_limit)){
#			if ((($taglist_line[1] <= $lower_limit) || (($taglist_line[1] >= $coverage_limits_lower[0])&&($taglist_line[1] <= $coverage_limits_lower[1]))) || ($taglist_line[2] <= $lower_limit)){
                                $ref = join ("\t", $taglist_line[0], $taglist_line[6]);
                                push (@read_ids_lower, $ref);
                        }elsif (($taglist_line[1] >= $upper_limit) || ($taglist_line[2] >= $upper_limit)){
                                $ref = join ("\t", $taglist_line[0], $taglist_line[6]);
                                push (@read_ids_upper, $ref);
                        }
                }
        }
        %ids = map { $_ => 1 } @read_ids_lower;
        my @unique_lower = keys %ids;
        %ids = map { $_ => 1 } @read_ids_upper;
        my @unique_upper = keys %ids;
        %ids = map { $_ => 1 } @read_ids_all;
        my @unique_all = keys %ids;
        for (@unique_all) {
                my @junk = split /\//;
                push (@reads, $junk[0]);
                @junk = split /\t/;
                push (@total_readlist, $junk[1]);
        }

        map { $count{$_}++ } @reads;
        map {if ($count{$_} == 2){ @readid_good = split /\t/; push(@readlist, "$readid_good[1]");} elsif ($count{$_} == 1) { push(@readid_bad, "$_");}} keys (%count);
        @reads = {};
        undef %count;

        for (@readlist){
                chomp;
                $current_id = $_;
                my @pairs_lower = grep { $_ =~ /$current_id/} @unique_lower;
                my @pairs_upper = grep{ $_ =~ /$current_id/} @unique_upper;
#                print "good id: $current_id\n";
                my $count_lower = scalar @pairs_lower;
                my $count_upper = scalar @pairs_upper;
#                print "count lower: $count_lower\n";
#                print "count upper: $count_upper\n";
                unless ((scalar @pairs_lower == 2) || (scalar @pairs_upper == 2)){
                        push (@readlist_good, "$current_id");
                }
        }
        for (@readid_bad){
                chomp;
                @unsorted = ();
                ($junk, $current_id) = split (/\t/);
                print "current ID: $current_id\n";
                @singleton = grep { $_ =~ /$current_id/} @total_readlist;
                for (@singleton){
                        chomp;
                        $tag = $_;
                        @taglist = grep { $_ =~ /$tag/} @readtaglist;
#                        print "taglist: @taglist\n";
                }
                for (@taglist) {
                        @taglist_line = split /\t/;
                        push(@unsorted, $taglist_line[1], $taglist_line[2]);
                        $max = max @unsorted;
                        $min = min @unsorted;
                }
#                print "unsorted: @unsorted\n";
		print "read mapping from $min to $max\n";
#                print "min: $min\n";
#                print "max: $max\n";

		if ($min <= $lower_limit){
                        print "orphan discarded! min<lowerlimit\n----------------------------------------------\n";
                }elsif ($max >= $upper_limit){
			print "orphan discarded! max>upperlimit\n----------------------------------------------\n";
		}elsif (($min >= $lower_main_limit) && ($max <= $upper_main_limit)){
			print "orphan discarded! lower_main_limit<min-max<upper_main_limit\n----------------------------------------------\n";
		}elsif (($min >= $elevated_cov_lower_start) && ($min <= $elevated_cov_lower_end - ($lower_limit / 2))){
			print "orphan discarded! increased_coverage_lower_start<min<increased_coverage_lower_end\n----------------------------------------------\n";
		}elsif (($max >= ($elevated_cov_upper_start + ($lower_limit / 2)))  && ($max <= $elevated_cov_upper_end)){
			print "orphan discarded! increased_coverage_upper_start<max<increased_coverage_upper_end\n----------------------------------------------\n";
		}else {
                        push(@singletons, "@singleton\n");
                        print "orphan resurrected! \n----------------------------------------------\n";
                }
#                print "contiglength: $contiglength\n";
        }

        for (@singletons){
                my @resurrection = split /\//;
                push (@readlist_good, $resurrection[0]);
        }
        for (@readlist_good){
                $_ =~ s/$/\/1\n/g;
                push(@readlist_proofread, $_);
                $_ =~ s/1$/2/g;
                push(@readlist_proofread, $_);
        }
        return @readlist_proofread;
}

sub standard_deviation {
        my(@numbers) = @_;
        #Prevent division by 0 error in case you get junk data
        return undef unless(scalar(@numbers));

        # Step 1, find the mean of the numbers
        my $total1 = 0;
       	foreach my $num (@numbers) {
		if (!$num){
			$num = 0;
		}
		$total1 += $num;
        }
        my $mean1 = $total1 / (scalar @numbers);
	push (my @stdev, "$mean1");
        # Step 2, find the mean of the squares of the differences between each number and the mean
        my $total2 = 0;
        foreach my $num (@numbers) {
		if (!$num){
			$num = 0;
		}
                $total2 += ($mean1-$num)**2;
        }
        my $mean2 = $total2 / (scalar @numbers);

        # Step 3, standard deviation is the square root of the above mean
        my $std_dev = sqrt($mean2);
	push (@stdev, "$std_dev");
        return @stdev;
}

sub assess_coverage{
        my $readtaglist_FH = $_[0];
        my @readtaglist;
        my $from =$_[1];
        my $to = $_[2];
	my $where = $_[3];
        my @taglist_line;
        my @coverage_array_lower;
        my @coverage_array_upper;
        my @read_ids_lower;
        my @read_ids_upper;
        my %ids;
        my @taglist;
        my @unsorted;
        my $min;
        my $max;
        my %coverage;
        my @allnums;
        my @coverage_change_position;
        my @coverage_limits;

#	print "assessing coverage from position $from to position $to\n";
	open (TAGLIST,"<$readtaglist_FH") or die $!;
	while (<TAGLIST>){
		push (@readtaglist, "$_");
	}
	close (TAGLIST);

        for (@readtaglist){
                @taglist_line = split /\t/;
                unless ($taglist_line[0] =~ /#/){
                        if ((($taglist_line[1] >= $from) && ($taglist_line[1] <= $to)) || (($taglist_line[2] >= $from) && ($taglist_line[2] <= $to))){
				push (@coverage_array_lower, "$_");
                                push (@read_ids_lower, $taglist_line[6]);
                        }
                }
        }
        %ids = map { $_ => 1 } @read_ids_lower;
        my @unique_lower = keys %ids;
        for (@unique_lower){
                my @current_id = $_;
                chomp;
                @unsorted = ();
                for (@current_id){
                        my $current_id = $_;
                        @taglist = grep { $_ =~ /$current_id/} @coverage_array_lower;
                }
                for (@taglist) {
                        @taglist_line = split /\t/;
                        push(@unsorted, $taglist_line[1], $taglist_line[2]);
                        $max = max @unsorted;
                        $min = min @unsorted;
                }
                my @nums = ($min .. $max);
                for (@nums){
                        push (@allnums, "$_");
                }
        }
	%coverage = map { $_ => 0 } @allnums;
        map { $coverage{$_}++ } @allnums;
        open (OUT,">out-$where.csv");

########## detecting coverage peak
        my $max_cov = 0;
	my $max_cov_position;
	my @cumulative_coverage;
	map { unless (!$coverage{$_}){print OUT "$_,$coverage{$_}\n"; push (@cumulative_coverage, "$coverage{$_}"); if ($coverage{$_} > $max_cov){$max_cov = $coverage{$_}; $max_cov_position = $_; }}} ($from..$to);
        my @average_coverage = &standard_deviation(@cumulative_coverage);
	my $coverage_factor = $max_cov / $average_coverage[0];
	open (OUT,">>out-$where.csv");
       	print OUT "\nmaximum coverage is $max_cov at position $max_cov_position\naverge coverage is: $average_coverage[0], sd: $average_coverage[1]\nfactor $coverage_factor\n";
	close (OUT);

######### detecting rapid changes in coverage
        for ($from..($to - 10)){
                my $position = $_;
                my $cov = $coverage{$position};
                unless (!$cov){
                        my @positions = ();
                        push (@positions, "$cov");
                        for (1 .. 10){
                                my $pos_plus = $position + $_;
                                if ($coverage{$pos_plus}){
                                        push (@positions, "$coverage{$pos_plus}");
                                }
                        }
			my @stdev = &standard_deviation(@positions);
                	if ($stdev[1] > 6.0){
                        	print "positions ($position): @positions -> stdev: $stdev[1]\n";
                        	push (@coverage_change_position, $position);
                	}elsif ($stdev[1] >= 4.5){
				print "positions ($position): @positions -> stdev: $stdev[1]\n";
                	}
                }
        }
        if (@coverage_change_position){
                print "positions with rapidly changing coverage detected: @coverage_change_position\n";
                my $start = min @coverage_change_position;
                my $end = max @coverage_change_position;
                push (@coverage_limits, "$start", "$end");
                print "set limits from $coverage_limits[0] to $coverage_limits[1]\n";
        }else{
                print "no irregularities in coverage detected\n";
                push (@coverage_limits, "0", "0");
        	return @coverage_limits;
        }


###### assessing whether coverage peak lies within putative conserved region, if yes accept prediction; if no, reject conserved region
	if (($coverage_factor >= 1.6) && (($coverage_limits[0] < $max_cov_position) && ( $max_cov_position < $coverage_limits[1]))){
		
		print "suspicious coverage peak detected within the predicted limits\n";
	}else {
		print "no coverage peak detected within predicted limits - rejecting limits\n";
		@coverage_limits = ();
		push (@coverage_limits, "0", "0");		
	}
	return @coverage_limits;	

}

sub check_ref_length{
        my $ref=$_[0];
        my $output_filename=$_[1];
        my $critical=$_[2];
	my $kmer = $_[3];
        my @header;
        my $header_count=0;
        my (@sequence,@temp_output,@final_output);
        my $full_sequence;
        open(REF,"<$ref") or die $!;
        while(<REF>){
                chomp;
                if ($_ =~ /^>/){
                        push(@header,$_);
#                       print "found header:\n$header[$header_count]\n";
                        $header_count++;
#                       print "header count: $header_count\n";
                        if (@sequence){
                                @temp_output=&finalize_sequence($critical,$header[-2],$kmer,@sequence);
                                for (@temp_output){
                                        push(@final_output,$_);
                                }
                        }
                        undef @sequence;
                }elsif ($_ =~ /[a-zA-Z]/){
#                       print "found sequence:\n$_\n";
                        push(@sequence,$_);
                }
        }
        @temp_output=&finalize_sequence($critical,$header[-1],$kmer,@sequence);
        for (@temp_output){
                push(@final_output,$_);
        }
#       print "result:\n";
        open (OUT,">$output_filename") or die $!;
        for(@final_output){
#               print "$_\n";
                print OUT "$_\n";
        }
        close REF;
        close OUT;

}

sub finalize_sequence{
        my $critical = shift;
        my $header = shift;
	my $kmer = shift;
        my $full_sequence=join("",@_);
        my $factor;
        my @output;
        if (!$critical){
                $factor=0;
        }else{
                $factor=ceil(length($full_sequence)/$critical);
        }
        if ($factor == 1){
                push(@output,$header);
                push(@output,$full_sequence);
        }else{ #too long
                print "\nreference is too long for mirabait to be handled in one go -> will be split into sub-sequences\n";
                $header=substr $header, 1;
                for (my $i=0; $i<$factor; $i++){
                        unless ((length(substr $full_sequence, $i*$critical, $critical+$kmer)-$kmer)<0){
                                push(@output,">sub$i\_" .$header);
                                push(@output,substr $full_sequence, $i*$critical, $critical+$kmer);
                        }
                }
        }
        return @output;
}

sub create_manifest {
	my ($iter, $sampleID, $refID, $mmode, $trim, $platform, $solexa_missmatches, $pair, $overhang, $reads, $ref, $redirect, $NFS_warn);
	$iter = $_[0];
	$sampleID = $_[1];
	$refID = $_[2];
	$mmode = $_[3];
	$trim = $_[4];
	$platform = $_[5];
	$solexa_missmatches = $_[6];
	$pair = $_[7];
	$overhang = $_[8];
	$reads = $_[9];
	$ref = $_[10];
	$redirect = $_[11];
	$NFS_warn = $_[12];

	if ($NFS_warn){
		$NFS_warn = ":cnfs=warn"
	}

	open (MANIFEST,">manifest.conf") or die $!;
	print MANIFEST "#manifest file for iteration $iter created by MITObim\n\nproject = $sampleID-$refID
	\njob = genome,$mmode,accurate
	\nparameters = -NW:mrnl=0$NFS_warn -AS:nop=1 $redirect $overhang $platform $trim -CO:msr=no $solexa_missmatches\n";
	my @technology = split("_",$platform);
	#-notraceinfo -
	if ($mmode eq "mapping"){
		print MANIFEST "\nreadgroup\nis_reference\ndata = $ref\nstrain = $refID\n";
	}
	print MANIFEST "\nreadgroup = reads\ndata = $reads\ntechnology = $technology[0]";
	if ($pair){
#		print MANIFEST "\nsegmentplacement = ---> <--- infoonly";
		print MANIFEST "\nsegmentplacement = ---> <--- exclusion_criterion";
		print MANIFEST "\ntemplatesize = -1 600 exclusion_criterion";
	}
	print MANIFEST "\nstrain = $sampleID\n";
	close MANIFEST;
}
sub clean {
        my $interval = shift;
        my $cur = shift;
        my $dir = $cur-$interval;
        my $path=abs_path;
        if (-d "$path/../iteration$dir"){
                print "\nnow removing directory iteration$dir\n";
                rmtree ("$path/../iteration$dir") or die $!;
        }
}

sub remove_unmapped_contigs{
	my $file = shift;
	my $out = shift;
	my ($head, $seq,$length);
	my @split;
	my ($dropped,$split,$trim)=(0,0,0);
	open (FASTA,"<$out") or die $!."\ncould not open $file for reading\n";
	while (<FASTA>){
		chomp;
		if ($_ =~ /^>/){
			$head=$_;
		}else {
			if ($_ =~ /X/){
				$dropped++;
				print "drop:\n$head\n$_\n";
			}else{
				$length=length($_);
				$_ =~ s/^N+//g;	#remove all Ns from the beginning of the sequence
				$_ =~ s/N+$//g;	#remove all Ns from the end of the sequence
				if ($length != length($_)){
					$trim++;
				}
				@split=split(/N{3,}/);	#split at every occurance of 3 or more Ns
			if (@split > 1){
					print "contig $head has been split into ".scalar @split." sub-contigs\n";
					$split++;
					for (my $i=0; $i<@split; $i++){
#						print "$head\_$i\n$split[$i]\n";
						$seq.="$head\_$i\n$split[$i]\n";
					}
				}else{
#					print "$head\n$_\n";
					$seq.="$head\n$_\n";
				}
			}
		}
	}		
	close FASTA;
	if (($dropped) || ($split) || ($trim)){
		print "creating backup of $out before changing it\n";
		copy $out,$file;
		open (OUT,">$out") or die $!."\ncould not open $out for writing\n";
		print OUT "$seq\n";
		close OUT;
	}

	if ($dropped){
		print "\nremoved $dropped contigs from baitfile because they did not receive any hits\n";
	}else{
		print "\nno need to remove any sequences from the baitfile\n";
	}

	if (!$split){
		print "\nno need to split any contigs\n";
	}
	if ($trim){
		print "\n$trim contigs had N's trimmed off their ends\n";
	}
}
