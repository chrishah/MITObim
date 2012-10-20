#! /usr/bin/perl
#
# MITObim - mitochondrial baiting and iterative mapping
# wrapper script version 1.1
# Author: Christoph Hahn, Oct. 2012
#
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Copy;

my $startiteration = 1;
my $enditeration = 1;
my $quick = 0;
my $noshow = 0;
my $help = 0;
my $strainname = 0;
my $paired = 0;
my $mode = 0;
my $refname = 0;
my $readpool = 0;
my $maf = 0;
my $optionstest = 0;

GetOptions (	"start=i" => \$startiteration,
		"end=i" => \$enditeration,
		"quick=s" => \$quick,
		"noshow!" => \$noshow,
		"strainname=s" => \$strainname,
		"paired" => \$paired,
		"denovo" => \$mode,
		"ref=s" => \$refname,
		"readpool=s" => \$readpool,
		"help!" => \$help,
		"maf=s" => \$maf) or die "Incorrect usage!\n";



print "\nMITObim - mitochondrial baiting and iterative mapping\n";
print "version 1.1\n";
print "author: Christoph Hahn, (c) 2012\n\n";

$optionstest = 1 if $help;
$optionstest = 1 if ($startiteration > $enditeration);
$optionstest = 1 if !$readpool;
unless ($quick){
	$optionstest = 1 if !$maf;
}
$optionstest = 1 if !$refname;

if ($optionstest ==1){
	print	"\nusage: ./MITObim.pl <parameters>\n";
	print 	"\nparameters:\n";
	print   "\n	-start <int>		iteration to start with, default=1\n";
	print   "	-end <int>		iteration to end with, default=1\n";
	print   "	-strain	<string>	strainname as used in initial MIRA assembly\n";
	print   "	-ref <string>		referencename as used in initial MIRA assembly\n";
	print   "	-readpool <PATH>	path to readpool in fastq format\n";
	print   "	-maf <PATH>		path to maf file from previous MIRA assembly\n";	
	print   "\noptional:\n";
	print   "\n	--denovo		runs MIRA in denovo mode, default: mapping\n";
	print   "	--pair			finds pairs after baiting, default: no\n";
	print	"	--quick <PATH>		starts process with initial baiting using provided fasta reference\n";
	print   "	--noshow		do not show output of MIRA modules\n";
	print   "	--help			shows this information\n";
	print	"\nexamples:\n";
	print   "./MITObim.pl -start 1 -end 5 -strain StrainX -ref Gthymalli-mt -readpool /PATH/TO/readpool.fastq -maf /PATH/TO/assembly.maf\n";
	print   "./MITObim.pl --quick /PATH/TO/reference.fasta -strain StrainY -ref Gthymalli-mt -readpool /PATH/TO/readpool.fastq\n";
	exit;
}
unless (((-e $maf)||($quick)) && (-e $readpool)){
        print "\nAre readpool AND maf files there?\n";
        exit;
}
if ($quick){
        unless (-e $quick){
		print "quick option selected but is the file there?\n";
		exit;
	}
	print "quick option selected! -maf option will be ignored (if given)\n";
	$maf = 0;
	$startiteration = 0;
}
if (!$mode){
	$mode = "mapping";
}else {
	$mode = "denovo";
}

print "\nAll paramters seem to make sense:\n";
print "startiteration: $startiteration\n";
print "enditeration: $enditeration\n";
print "strainname: $strainname\n";
print "refname: $refname\n";
print "readpool: $readpool\n";
print "maf: $maf\n";
print "quick: $quick\n";
print "paired: $paired\n";
print "denovo: $mode\n";
print "noshow: $noshow\n";

print "\nStart MITObim\n";

my @iteration = ($startiteration .. $enditeration);
foreach (@iteration){
	chomp;
	mkdir "iteration$_" or die "$!";
	chdir "iteration$_" or die "$!";
	print "\n============\n";
	print " ITERATION$_\n";
	print "============\n\n";
	if ($maf){
		print "\nrecover backbone by running convert_project on maf file\n";

		my @convert_project_output = (`convert_project -f maf -t fasta -A "SOLEXA_SETTINGS -CO:fnicpst=yes" $maf tmp 2>&1`);
		my $convert_exit = $? >> 8;
		unless ($noshow){
			print "@convert_project_output\n";
		}
		unless ($convert_exit == 0){
			print "\nconvert_project seems to have failed - see detailed output above\n";
			exit;
		}
		open(FH1,"<tmp_$strainname.unpadded.fasta") or die "$!";
		open(FH2,">$strainname-$refname\_backbone_in.fasta") or die "$!";
		while (<FH1>) {
			$_ =~ s/@/N/g;
			print FH2 $_; 
		}
		close(FH1);
		close(FH2);
		unlink glob ("tmp*");
	}
	MIRABAIT:
	unless ($maf){
		print "\nquick option baits reads from provided reference in iteration0\n";
		copy("$quick", "$strainname-$refname\_backbone_in.fasta") or die "copy failed: $!";		
	}
	print "\nfishing readpool using mirabait\n";
	
	my @mirabait_output = (`mirabait -k 31 -n 1 $strainname-$refname\_backbone_in.fasta $readpool $strainname-$refname\_in.solexa 2>&1`);
	my $mirabait_exit = $? >> 8;
	unless ($noshow){
		print "@mirabait_output\n";
	}
	unless ($mirabait_exit == 0){
	        print "\nmirabait seems to have failed - see detailed output above\n";
	        exit;
	}
	
	FINDPAIRS:
	
	unless (!$paired){
		print "\nfind pairs to baited reads\n";
		open(FH1,"<$strainname-$refname\_in.solexa.fastq");
		open(FH2,">list");
		my %hash;
		my @temp;
		my $key;
		my $val;

		while (<FH1>) {
			if ($_ =~ /^\@HWI/){
                		$_ =~ s/@//g;
				($key, $val) = split /\//;
		        	$hash{$key} .= exists $hash{$key} ? ",$val" : $val;
			}
		}
		for (keys %hash){
			$_ =~ s/$/\/1/g;
			print FH2 "$_\n";
			$_ =~ s/1$/2/g;
	        	print FH2 "$_\n";
		}
		close(FH1);
		close(FH2);
		my @convert_project_output = (`convert_project -f fastq -t fastq -n list $readpool $strainname-$refname\_in.solexa 2>&1`);
		my $convert_exit = $? >> 8;
		unless ($noshow){
			print "@convert_project_output\n";
		}
		unless ($convert_exit == 0){
	        	print "\nconvert_project seems to have failed - see detailed output above\n";
	        	exit;
		}
	}
	unlink("list");
	MIRA:
	print "\nrunning assembly using MIRA v3.4\n\n";
	my @MIRA_output = (`mira --project=$strainname-$refname --job=$mode,genome,accurate,solexa "--noclipping -CL:pec=no" -MI:somrnl=0 -AS:nop=1 -SB:bsn=$refname:bft=fasta:bbq=30 SOLEXA_SETTINGS -CO:msr=no -GE:uti=$paired -SB:dsn=$strainname 2>&1`);
	my $MIRA_exit = $? >> 8;
	unless ($noshow){
		print "@MIRA_output\n";
	}
	unless ($MIRA_exit == 0){
		print "\nMIRA seems to have failed - see detailed output above\n";
		exit;
	}
	my @path = abs_path;
	push (@path, "/$strainname-$refname\_assembly/$strainname-$refname\_d_results/$strainname-$refname\_out.maf");
	$maf = join("",@path);
	unless (-e $maf){
		print "maf file is not there \n";
		exit;
	}
	chdir ".." or die "Failed to go to parent directory: $!";
}
print "\nsuccessfully completed $enditeration iterations with MITObim!\n\n";
