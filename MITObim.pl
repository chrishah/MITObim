#! /usr/bin/perl
#
# MITObim - mitochondrial baiting and iterative mapping
# wrapper script version 1.3
# Author: Christoph Hahn, Oct. 2012
#
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Copy;
use List::Util qw< min max >;
use POSIX qw(strftime);

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
my $proofreading = 0;
my $readlength = 0;
my $insertsize = 0;
#my $optionstest = 0;
my $miramode;
my @reads = ();
my %hash;
my $key;
my $val;
my @output = ();
my $exit = ();
my @path = ();
my $shme = "";
my $MM = 0;
my $current_contiglength;
my @contiglengths;
my $current_number_of_reads;
my @number_of_reads;
my $USAGE = 	"\nusage: ./MITObim.pl <parameters>\n
	 	\nparameters:\n
		-start <int>		iteration to start with, default=1
		-end <int>		iteration to end with, default=1
		-strain	<string>	strainname as used in initial MIRA assembly
		-ref <string>		referencename as used in initial MIRA assembly
		-readpool <PATH>	path to readpool in fastq format
		-maf <PATH>		path to maf file from previous MIRA assembly\n
		\noptional:\n
		--denovo		runs MIRA in denovo mode, default: mapping
		--pair			finds pairs after baiting, default: no
		--quick <PATH>		starts process with initial baiting using provided fasta reference
		--noshow		do not show output of MIRA modules
		--help			shows this information
		--proofread		applies proofreading
		--readlength <int>	read length of illumina library, default=150
		--insert <int>		insert size of illumina library, default=300
		\nexamples:\n
		./MITObim.pl -start 1 -end 5 -strain StrainX -ref Gthymalli-mt -readpool /PATH/TO/readpool.fastq -maf /PATH/TO/assembly.maf
		./MITObim.pl --quick /PATH/TO/reference.fasta -strain StrainY -ref Gthymalli-mt -readpool /PATH/TO/readpool.fastq\n";

my $PROGRAM = "\nMITObim - mitochondrial baiting and iterative mapping\n";
my $VERSION = "version 1.9\n";
my $AUTHOR = "author: Christoph Hahn, (c) 2012\n\n";

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
		"maf=s" => \$maf,
		"proofreading!" => \$proofreading,
		"readlength=i" => \$readlength,
		"insertsize=i" => \$insertsize) or die "Incorrect usage!\n$USAGE";



print $PROGRAM; 
print $VERSION; 
print $AUTHOR; 

print $USAGE and exit if $help;
print $USAGE and exit if ($startiteration > $enditeration);
print $USAGE and exit if !$readpool;
unless ($quick){
        print $USAGE and exit if !$maf;
}
print $USAGE and exit if !$refname;

unless (((-e $maf)||($quick)) && (-e $readpool)){
        print "\nAre readpool AND maf files there?\n";
        exit;
}
if (!$readlength){
	$readlength = 150;
}
if (!$insertsize){
	$insertsize = 300;
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
	$miramode = "mapping";
}else {
	$miramode = "denovo";
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
print "denovo: $mode (mapping=0, denovo=1)\n";
print "noshow: $noshow\n";
print "proofread: $proofreading\n";
print "readlength: $readlength\n";
print "insertsize: $insertsize\n";
if ($proofreading){
	print "number of allowed missmatches in proofreading assembly: $MM\n";
}

print "\nStarting MITObim \n";

my @iteration = ($startiteration .. $enditeration);
foreach (@iteration){
	chomp;
	my $currentiteration = $_;
	mkdir "iteration$currentiteration" or die $!;
	chdir "iteration$currentiteration" or die $!;
	print "\n==============\n";
	print " ITERATION $currentiteration\n";
	print "==============\n";
	print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
#	if (($proofreading) && ($currentiteration != 0)){
	if ($proofreading){
		$shme = "-AL:shme=$MM";
	} 
		
	if ($maf){
		print "\nrecover backbone by running convert_project on maf file\n";

		@output= (`convert_project -f maf -t fasta -A "SOLEXA_SETTINGS -CO:fnicpst=yes" $maf tmp 2>&1`);
		$exit = $? >> 8;
		unless ($noshow){
			print "@output\n";
		}
		unless ($exit == 0){
			print "\nconvert_project seems to have failed - see detailed output above\n";
			exit;
		}
		if ( ((($mode) && ($currentiteration > 1)) && (!$quick)) || ((($mode) && ($currentiteration >= 1)) && ($quick)) ){
			open(FH1,"<tmp_default.unpadded.fasta") or die "$!";
		}else{
			open(FH1,"<tmp_$strainname.unpadded.fasta") or die "$!";
		}
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
		print "\nquick option baits reads from provided reference in iteration 0\n";
		copy("$quick", "$strainname-$refname\_backbone_in.fasta") or die "copy failed: $!";		
	}
	print "\nfishing readpool using mirabait\n";
	
	@output = (`mirabait -k 31 -n 1 $strainname-$refname\_backbone_in.fasta $readpool $strainname-$refname\_in.solexa 2>&1`);
	$exit = $? >> 8;
	unless ($noshow){
		print "@output\n";
	}
	unless ($exit == 0){
	        print "\nmirabait seems to have failed - see detailed output above\n";
	        exit;
	}
	
	FINDPAIRS:
	
	unless (!$paired){
		print "\nfind pairs to baited reads\n";
		open(FH1,"<$strainname-$refname\_in.solexa.fastq") or die $!;
		open(FH2,">list");

		while (<FH1>) {
			if ((substr($_, 0, 1) == "@") && (substr($_, 1, 1) != "@")){
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
		@output = (`convert_project -f fastq -t fastq -n list $readpool $strainname-$refname\_in.solexa 2>&1`);
		$exit = $? >> 8;
		unless ($noshow){
			print "@output\n";
		}
		unless ($exit == 0){
	        	print "\nconvert_project seems to have failed - see detailed output above\n";
	        	exit;
		}
	}
	unlink("list");
	
	MIRA:
	print "\nrunning assembly using MIRA v3.4\n\n";
	@output = (`mira --project=$strainname-$refname --job=$miramode,genome,accurate,solexa "--noclipping -CL:pec=no" -MI:somrnl=0 -AS:nop=1 -SB:bsn=$refname:bft=fasta:bbq=30 SOLEXA_SETTINGS -CO:msr=no -GE:uti=$paired $shme -SB:dsn=$strainname 2>&1`);
	$exit = $? >> 8;
	unless ($noshow){
		print "@output\n";
	}
	unless ($exit == 0){
		print "\nMIRA seems to have failed - see detailed output above\n";
		exit;
	}

	@path = abs_path;
        push (@path, "/$strainname-$refname\_assembly/$strainname-$refname\_d_results/$strainname-$refname\_out.maf");
        $maf = join("",@path);
        unless (-e $maf){
                print "maf file is not there \n";
                exit;
        }
	$current_contiglength = &get_contig_length("$strainname-$refname\_assembly/$strainname-$refname\_d_info/$strainname-$refname\_info_contigstats.txt");
	$current_number_of_reads = &get_number_of_reads("$strainname-$refname\_assembly/$strainname-$refname\_d_info/$strainname-$refname\_info_contigstats.txt");

	PROOFREAD:
#	if (($proofreading) && ($currentiteration >= 1)){
	if ($proofreading){
		print "\n Proofreading\n\n";
		my $contigreadlist = "$strainname-$refname\_assembly/$strainname-$refname\_d_info/$strainname-$refname\_info_contigreadlist.txt";
		my $readtaglist = "$strainname-$refname\_assembly/$strainname-$refname\_d_info/$strainname-$refname\_info_readtaglist.txt";
		
#		print "assessing coverage in lower region from position 0 to 500\n";
		my @coverage_limits_lower = &assess_coverage($readtaglist, 0, (2*$insertsize), "lower");
#		print "assessing coverage in upper region from position " . ($current_contiglength-500) . " to $current_contiglength\n";
		my @coverage_limits_upper = &assess_coverage($readtaglist, ($current_contiglength - (2*$insertsize)), ($current_contiglength), "upper");
		
#		mkdir "tmp" or die $!;
#        	chdir "tmp" or die $!;
#		copy ("../$strainname-$refname\_backbone_in.fasta", "$strainname-$refname-1MM_backbone_in.fasta") or die "copy failed: $!";
#		copy ("../$strainname-$refname\_in.solexa.fastq", "$strainname-$refname-1MM_in.solexa.fastq") or die "copy failed: $!";
        	
#		print "\nrunning second instance of MIRA during proofreading\n\n";
#		@output = (`mira --project=$strainname-$refname-1MM --job=$miramode,genome,accurate,solexa "--noclipping -CL:pec=no" -MI:somrnl=0 -AS:nop=1 -SB:bsn=$refname:bft=fasta:bbq=30 SOLEXA_SETTINGS -CO:msr=no -GE:uti=yes:tismin=100:tismax=600 -AL:shme=3 -SB:dsn=$strainname 2>&1`);
#        	$exit = $? >> 8;
#        	unless ($noshow){
#                	print "@output\n";
#        	}
#        	unless ($exit == 0){
#                	print "\nMIRA seems to have failed - see detailed output above\n";
#                	exit;
#        	}
#		chdir ".." or die $!;
#		my $contigreadlist_1MM = "tmp/$strainname-$refname-1MM_assembly/$strainname-$refname-1MM_d_info/$strainname-$refname-1MM_info_contigreadlist.txt";
		print "\nScreening orphan reads and discarding potentially dubious reads\n";	
		open(OUT,">list");
#		print OUT &proofread($contigreadlist, $contigreadlist_1MM);
		print OUT &proofread($contigreadlist, $readtaglist, $current_contiglength, $coverage_limits_lower[0], $coverage_limits_lower[1], $coverage_limits_upper[0], $coverage_limits_upper[1], (1.3*$readlength), (2*$insertsize));	
		close(OUT);	

		print "\ngenerating proofread readpool\n";
		@output = (`convert_project -f fastq -t fastq -n list $strainname-$refname\_in.solexa.fastq $strainname-$refname-proofread\_in.solexa 2>&1`);
		$exit = $? >> 8;
		unless ($noshow){
			print "@output\n";
		}
		unless ($exit == 0){
			print "\nconvert_project seems to have failed - see detailed output above\n";
			exit;
		}	
		copy("$strainname-$refname\_backbone_in.fasta", "$strainname-$refname-proofread\_backbone_in.fasta") or die "copy failed: $!";	

		print "\nrunning proofread assembly using MIRA v3.4\n\n";
        	@output = (`mira --project=$strainname-$refname-proofread --job=$miramode,genome,accurate,solexa "--noclipping -CL:pec=no" -MI:somrnl=0 -AS:nop=1 -SB:bsn=$refname:bft=fasta:bbq=30 SOLEXA_SETTINGS -CO:msr=no -GE:uti=yes:tismin=100:tismax=600 $shme -SB:dsn=$strainname 2>&1`);
        	$exit = $? >> 8;
        	unless ($noshow){
                	print "@output\n";
        	}
        	unless ($exit == 0){
                	print "\nMIRA seems to have failed - see detailed output above\n";
                	exit;
        	}
		
		@path = abs_path;
		push (@path, "/$strainname-$refname-proofread\_assembly/$strainname-$refname-proofread\_d_results/$strainname-$refname-proofread\_out.maf");
		$maf = join("",@path);
		unless (-e $maf){
			print "maf file is not there \n";
			exit;
		}
		$current_contiglength = &get_contig_length("$strainname-$refname-proofread_assembly/$strainname-$refname-proofread_d_info/$strainname-$refname-proofread_info_contigstats.txt");
		$current_number_of_reads = &get_number_of_reads("$strainname-$refname-proofread_assembly/$strainname-$refname-proofread_d_info/$strainname-$refname-proofread_info_contigstats.txt");

	}
	push (@number_of_reads, "$current_number_of_reads");
	push (@contiglengths, "$current_contiglength");
	print "readpool contains $current_number_of_reads reads\n";
	print "contig length: $contiglengths[-1]\n";
	if ($number_of_reads[-2]){
		if ($number_of_reads[-2] == $number_of_reads[-1]){
			print "\nMITObim has reached a stationary read number after $currentiteration iterations!!\n";
			print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
			exit;
		}
	}
#	while ($proofreading) {
#		if ($contiglengths[-2]){
#			if ($contiglengths[-1] == $contiglengths[-2]){
#				print "\nMITObim has reached a stationary contig length after $currentiteration iterations!!\n";
#				print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
#				exit;
#			}
#		}
#	}
	chdir ".." or die "Failed to go to parent directory: $!";
}
print "\nsuccessfully completed $enditeration iterations with MITObim! " . strftime("%b %e %H:%M:%S", localtime) . "\n\n";

#
#
###SUBROUTINES
#
#
#

sub get_contig_length{
	my $contig = $_[0];
	my @contiglength;
	open (CONTIGSTATS,"<$contig") or die $!;
	while (<CONTIGSTATS>){
		unless ($_ =~ /#/){
			@contiglength = split /\t/;
		} 
	}
	close (CONTIGSTATS);		
	return $contiglength[1];
}

sub get_number_of_reads{
	my $contig = $_[0];
	my @contigstats;
	open (CONTIGSTATS,"<$contig") or die $!;
	while (<CONTIGSTATS>){
		unless ($_ =~ /#/){
			@contigstats = split /\t/;
		} 
	}
	close (CONTIGSTATS);		
	return $contigstats[3];
}

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
#                if (($min <= $lower_limit) || ($max >= $upper_limit) || ($min >= $lower_main_limit) || ($max <= $upper_main_limit)) {

		if ($min <= $lower_limit){
                        print "orphan discarded! min<lowerlimit\n----------------------------------------------\n";
                }elsif ($max >= $upper_limit){
			print "orphan discarded! max>upperlimit\n----------------------------------------------\n";
		}elsif (($min >= $lower_main_limit) && ($max <= $upper_main_limit)){
			print "orphan discarded! lower_main_limit<min-max<upper_main_limit\n----------------------------------------------\n";
#		}elsif ((($min >= $elevated_cov_lower_start) && ($min <= $elevated_cov_lower_end)) || (($max >= $elevated_cov_lower_start) && ($max <= $elevated_cov_lower_end)) || (($min < $elevated_cov_lower_start) && ($max > $elevated_cov_lower_end))){
		}elsif (($min >= $elevated_cov_lower_start) && ($min <= $elevated_cov_lower_end - ($lower_limit / 2))){
#		}elsif ( (($max >= $elevated_cov_lower_start) && ($max <= $elevated_cov_lower_end)) || (($min < $elevated_cov_lower_start) && ($max > $elevated_cov_lower_end))){
			print "orphan discarded! increased_coverage_lower_start<min<increased_coverage_lower_end\n----------------------------------------------\n";
#		}elsif ((($max >= $elevated_cov_upper_start) && ($max <= $elevated_cov_upper_end)) || (($min >= $elevated_cov_upper_start) && ($min <= $elevated_cov_upper_end)) || (($min < $elevated_cov_upper_start) && ($max > $elevated_cov_upper_end))){
		}elsif (($max >= ($elevated_cov_upper_start + ($lower_limit / 2)))  && ($max <= $elevated_cov_upper_end)){
#		}elsif ( (($min >= $elevated_cov_upper_start) && ($min <= $elevated_cov_upper_end)) || (($min < $elevated_cov_upper_start) && ($max > $elevated_cov_upper_end))){
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
#	return $std_dev;
}

sub assess_coverage{
        my $readtaglist_FH = $_[0];
#        my $contiglength = $_[1];
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



#		my $max_stdev = 0;
#		my $max_stdev_position = 0;
#		print "\nunusually high coverage detected in $where region! coverage peaks at position $max_cov_position\n";
#		for ($from..($max_cov_position - 5)){
##		for ($from..($to - 5)){
#                	my $position = $_;
#                	my $cov = $coverage{$position};
#                	unless (!$cov){
#                        	my @positions = ();
#                        	push (@positions, "$cov");
#                        	for (1..4){
#                                	my $pos_plus = $position + $_;
#                                	if ($coverage{$pos_plus}){
#						push (@positions, "$coverage{$pos_plus}");
#                        		}
#				}
#                		my @stdev = &standard_deviation(@positions);
#				if ($stdev[1] > $max_stdev){
#					$max_stdev = $stdev[1];
#					$max_stdev_position = $position;
#				}
#				if ($stdev[1] >= 2.0){
#					print "positions ($position): @positions -> stdev: $stdev[1]\n";
#				}
#
#                	if ($stdev[1] > 2.5){
#                        	print "positions ($position): @positions -> stdev: $stdev[1]\n";
#                        	push (@coverage_change_position, $position);
#			}elsif ($stdev[1] >= 2.0){
#                       	 print "positions ($position): @positions -> stdev: $stdev[1]\n";
#                	}
#                	}
#        	}
#		push (@coverage_limits, "$max_stdev_position");
#		$max_stdev = 0;
#                $max_stdev_position = 0;
#		for (reverse ($max_cov_position..($to - 5))){
#			my $position = $_;
#			my $cov = $coverage{$position};
#			unless (!$cov){
#                                my @positions = ();
#                                push (@positions, "$cov");
#                                for (1..4){
#                                        my $pos_plus = $position + $_;
#                                        if ($coverage{$pos_plus}){
#                                                push (@positions, "$coverage{$pos_plus}");
#                                        }
#                                }
#                                my @stdev = &standard_deviation(@positions);
#                                if ($stdev[1] > $max_stdev){
#                                        $max_stdev = $stdev[1];
#                                        $max_stdev_position = $position;
#                                }
#				if ($stdev[1] >= 2.0){
#					print "positions ($position): @positions -> stdev: $stdev[1]\n";
#				}
#			}
#		}
#		push (@coverage_limits, "$max_stdev_position");
#		print "conserved region predicted from position $coverage_limits[0] to $coverage_limits[1]\n";
#	}else {
#		print "\nno indication of increased coverage detected in $where region!\n\n";
#		push (@coverage_limits, "0", "0");
#	}
#        if (@coverage_change_position){
#                print "positions with rapidly changing coverage detected: @coverage_change_position\n";
#                my $start = min @coverage_change_position;
#                my $end = max @coverage_change_position;
#                push (@coverage_limits, "$start", "$end");
#                print "set limits from $coverage_limits[0] to $coverage_limits[1]\n";
#        }else{
#                print "no irregularities in coverage detected\n";
#        	push (@coverage_limits, "0", "0");
#	}
#        return @coverage_limits;
}

