#! /usr/bin/perl
#
# script to extract per base coverage info for unpadded
# assembly result from MIRA wiggle file
# part of the MITObim pipeline
# c.hahn@hull.ac.uk
# November 2015

use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);

my ($prefix, $wiggle, $window) = ("", "", 1);
my @files = ();
my %wiggle = ();
my %padded = ();
my %unpadded = ();
my ($instep, $inspan) = (0,0);
my $help;

my $command = $0;
for (@ARGV){
        $command .= " $_";
}

my $USAGE="\n	USAGE: get_wiggle.pl <options>

		The default wiggle file produced by MIRA matches with the *padded* result file(s). This produces a wiggle file that matches the desired *unpadded* fasta file.

		options:
		-p, --prefix	<string>	Prefix for MIRA assembly files (mandatory), expects to find prefix.unpadded.fasta and prefix.padded.fasta
		-w, --wiggle	<string>	Wiggle file (optional), if not specified script will expect to find a file *.wig in the same directory as assembly results
		-h, --help			Show this help text.\n\n";

GetOptions ( 	"prefix=s" => \$prefix,
                "wiggle=s" => \$wiggle,
		"help!" => \$help,
		"step=i" => \$window) or die "Incorrect usage!\n$USAGE";

print $USAGE and exit if $help;
print $USAGE and exit if !$prefix;

@files = &check_wigfiles($prefix, $wiggle);
%wiggle = &parse_wig($files[2]);
%padded = &parse_fasta($files[0]);
%unpadded = &parse_fasta($files[1]);
($instep, $inspan) = &find_step_span(\%wiggle);
&adjust_wiggle($inspan, $instep, \%wiggle, \%padded);
for my $key (keys %padded){
	print "padded - $key: ".(length($padded{$key}))."\n";
	print "unpadded - $key: ".(length($unpadded{$key}))."\n";
	print "wiggle - $key: ".(scalar(@{$wiggle{$key}})-1)."\n";

	&output_wiggle($prefix, @{$wiggle{$key}});
}

sub output_wiggle {
	my $pre = shift;
	my @wig = @_;
	my $wig_head=shift(@wig);
	my @wig_heads=split(/\|/, $wig_head);
	open(my $OUT_FH, '>', $prefix.".unpadded.wig") or die $!;
        print $OUT_FH "$wig_heads[0]\n$wig_heads[1]\n";
        for (@wig){
                print $OUT_FH $_."\n";
        }
	close($OUT_FH);
}

sub adjust_wiggle {
	my $inwindow = shift;
	my $instep = shift;
	my ($wig_ref, $pad_ref) = @_;
	my $wig_index = 1;
	my @new_wig = ();
	foreach (keys %$wig_ref){
#		print "$_\n";
#		print "@{${$wig_ref}{$_}}\n";
#		print "${$pad_ref}{$_}\n";
		for (my $i=0; $i<length(${$pad_ref}{$_}); $i++){
			if ($i % $instep == 0){
				my $current = substr(${$pad_ref}{$_}, $i, $inwindow);
				$current =~ s/\*//g;
				my $len = length($current);
				my $cov = ${${$wig_ref}{$_}}[$wig_index];
#				print "$i\t$current\t$cov\t$len\t".($cov x $len)."\n";
				push(@new_wig, ("$cov") x $len); #$cov for (1..$len);
				$wig_index++;
			}
		}
		print "\n".scalar(@new_wig)."\n";
	}
}

sub find_step_span {
	my ($wig_ref) = @_;
	my ($span, $step) = (0,0);
        for my $key (keys %$wig_ref){
		my @temp = split(" ", ${${$wig_ref}{$key}}[0]);
		for (@temp){
			if ($_ =~ /^span=/){
				$span = $_;
				$span =~ s/span=//;
			}elsif ($_ =~ /step=/) {
				$step = $_;
				$step =~ s/step=//;
			}
		}
		last;		
	}
	return ($step, $span);
}
sub parse_fasta {
	my $padded = shift;
	my %seq_hash = ();
	my ($head, $seq) = ("","");
	my @valid = ();
	open(FH, $padded) or die "$!";
	for (<FH>){
		chomp;
		if ($_ =~ /^>/){
			$_ =~ s/^>//;
			$head = $_;
			$seq_hash{$head} = "";
		}else{
			$seq_hash{$head} .= $_;
		}
#		print "$seq_hash{$_}\n";
	}
	return %seq_hash;
}

sub parse_wig {
	my $wiggle_file = shift;
	my $head = "";
	my $key = "";
	my %wiggle = ();
	open(FH, $wiggle_file) or die "$!";
	for (<FH>){
		chomp;
		unless ($_ =~ /^[0-9]+\z/){
#			print "$_\n";
			$head.=$_."|";
#			print "\n$head\n";
			$key = "";
		}else{
			if (!$key){
				my @temp = split(" ", $head);
					for (@temp){
						if ($_ =~ /^chrom=/){
							$_ =~ s/chrom=//;
							$key = $_;
#							print $key."\n";
						}
					}
#				$key = $head;
				$head = substr($head,0,-1);#~ s/|$//g;
				$wiggle{$key} = [$head];
				$head = "";
			}
#			print "$_\n";
			push(@{$wiggle{$key}}, $_);
#			print scalar(@{$wiggle{$key}})."\n";
		}
	}
	return %wiggle;
}

sub check_wigfiles{
	my $pre = shift;
	my $wig = shift;
	my $path = "";
	my @wig = ();
	my @files = ();
	if ((-e $pre.".unpadded.fasta") and (-e $pre.".padded.fasta")){
		push(@files, abs_path($pre.".padded.fasta"), abs_path($pre.".unpadded.fasta"));
	}else {
		print "can't find any sequence file matching your prefix\n";
		exit;
	}

	if ($wiggle){
		if (-e $wiggle){
#			print "ok - found wiggle file\n";
			push(@files, abs_path($wiggle));
		}else {
			print "specified wiggle file is not present\n";
			exit;
		}
	} else {
#		print "Trying to find wiggle file\n";
		$path = abs_path($pre.".unpadded.fasta");
		$path =~ s/\/${pre}.unpadded.fasta//;
		@wig = glob("$path/*.wig");
		unless (scalar(@wig)==1){
			print "either no or multiple wig files present - abort\n";
			exit;
		}else{
			push(@files, $wig[0]);
		}
	}
	
	print "All necessary files found:\n";
	for (@files){
		print "$_\n";
	}

	return @files;
}
#less Tthy-Salpinus_out_Tthy.padded.fasta | grep ">" -v | perl -ne 'chomp; @a=split(""); for (@a){print "$_\n"};' | perl -ne 'print "$_"; if ($_ =~ /\*/){ $count++}; if (($. % 4 == 0)&&$count){print "count is $count\n"; $count=0}' | grep "count is" |les
