#!/usr/bin/python

import sys,warnings
import argparse

VERSION="0.4"
DESCRIPTION='''
circules - checks for circularity in nucleotide sequences - version: v.%s

# Disclaimer: #
The script is currently in the beta phase. Any feedback is much appreciated.''' %VERSION
kmers = []
dic = {}
seqs = {}
auto=True

parser = argparse.ArgumentParser(description=DESCRIPTION, prog='circules.py', 
	formatter_class=argparse.RawDescriptionHelpFormatter, 
	epilog='''   examples: 
	# check for circularity using a k-mer length of 31 - returns suggestion(s) for clipping points in putative circular sequences.
	circules.py -f test.fasta -k 31
	
	# check for circularity using k-mers 31-37, stepsize 2bp - returns suggestion(s) for clipping points in putative circular sequences.
	circules.py -f test.fasta -k 31-37 -s 2
	
	# check for circularity using k-mers 31-37, stepsize 2bp. If length of suggested circular sequence is within +-10 percent of the expected length (15000 bp)
	# the clipped sequence will be written to a file called 'mytest.circular.fasta'. In addition a file 'mytest.for-testing.fasta' will be written. It contains
	# a 1000 bp sequence obtained by joining the first/last 500 bp of the proposed clipped circular sequence for additional evaluation if required.
	circules.py -f test.fasta -k 31-37 -s 2 -e 15000 -p mytest
	''')

parser.add_argument("-f", "--fasta", help="fasta file containing the sequence to be evaluate.", metavar="<FILE>", action="store")
parser.add_argument("-p", "--prefix", help="prefix for output files (default = 'circular').", type=str, metavar="<STR>", action="store", default='circular')
parser.add_argument("-e", "--exp_length", help="expected length (in bp) of circular molecule. If a candidate of length expected+-10%% is found, sequence will be clipped and written to file 'prefix.circular.fasta'.", type=int,  metavar="<INT>", default='0', action="store")
parser.add_argument("-t", "--length_tolerance", help="length tolerance. Candidate fragments must have a length of 'expected length +/- t * expected length'. Default = 0.1.", type=float,  metavar="<FLOAT>", default='0.1', action="store")
parser.add_argument("-k", "--kmer", help="kmer size. single number (default = 31) or range (e.g. 31-35).", type=str, default='31', action="store")
parser.add_argument("-s", "--kmer_step", help="kmer step size (default = 2).", type=int, default='2', action="store")
#parser.add_argument("-r", "--readpool", help="path to fastq reads to be used for evaluating circularity.", type=str, metavar="<FILE>", action="store")
parser.add_argument("-v", "--verbose", help="turn verbose output on.", action="store_true")
parser.add_argument("--version", action="version", version=VERSION) 
args = parser.parse_args()

if len(sys.argv) < 2:	#if the script is called without any arguments display the usage
    print
    parser.print_usage()
    print
    sys.exit(1)

if '-' in args.kmer:
	kmers = args.kmer.strip().split("-")
	kmers.append(args.kmer_step)
	for i in range(len(kmers)):
		kmers[i] = int(kmers[i])
else:
	kmers = [int(args.kmer),int(args.kmer)+1, 1]


fa = open(args.fasta,'r')
for l in fa:
	if l.startswith('>'):
		current = l.strip().replace('>','')
		seqs[current] = ''
	else:
		seqs[current]+=l.strip()

if len(seqs) > 1:
	print "\nThe file '%s' contains %i sequences. The circules script is currently designed to only handle a single fasta sequence at a time - Sorry!\n" %(args.fasta, len(seqs))
	sys.exit()

for k in range(kmers[0],kmers[1],kmers[2]):
	clips = {}
	print "\nCollecting %s-mers .." %k,
	for s in seqs:
		i=0
		while i <= (len(seqs[s])-k):
			if not seqs[s][i:i+k][-1] in 'acgtACGT':
				i+=k
				continue
			if not str(seqs[s][i:i+k]) in dic:
				dic[str(seqs[s][i:i+k])] = []
			dic[str(seqs[s][i:i+k])].append(i)
			i+=1
		print "Done!"
	print "Evaluation %s-mers .." %k,
	for s in seqs:
		i=0
		while i <= (len(seqs[s])-k):
			if not seqs[s][i:i+k][-1] in 'acgtACGT':
				i+=k
				continue
			if len(dic[str(seqs[s][i:i+k])]) == 2:
				length = dic[str(seqs[s][i:i+k])][1]-dic[str(seqs[s][i:i+k])][0]
				if not length in clips:
					clips[length] = []
				clips[length].append(dic[str(seqs[s][i:i+k])])
			i+=1
		print "Done!\n"
	suggest = []
	if len(clips) > 1:
		print "Found %i possibilites for circularity:" %len(clips)
		for l in sorted(clips):
			print "\t-> suggested circular length: %i (supported by %i %i-mers)" %(l, len(clips[l])/2, k)
#			if l < float(args.exp_length)*1.1 and l > float(args.exp_length)*0.9:
#				suggest.append(l)

	print ""

if not clips:
	print "Did not find any candidates for circularity\n"
	sys.exit()

if args.exp_length:
	print "Evaluating candidate lengths - specified expected range is: %s - %s" %(str(float(args.exp_length)*(1+args.length_tolerance)), str(float(args.exp_length)*(1-args.length_tolerance)))
	for l in sorted(clips):
		if l < float(args.exp_length)*(1+args.length_tolerance) and l > float(args.exp_length)*(1-args.length_tolerance):
			pass
			#print "%i is within +-10%% of specified expected length: %i" %(l, args.exp_length)
		else:
			del(clips[l])

	if clips:
		print "\nFound candidates within the specified length range"
		for s in seqs:
			for l in sorted(clips):
				#writing out clipped sequence
				print "length: %s - clip points: %i - %i -> writing sequence to '%s.circular.%s.fasta'\n" %(str(l), int(clips[l][0][0]), int(clips[l][0][1]), args.prefix, str(l))
				sequence = seqs[s][int(clips[l][0][0]):int(clips[l][0][1])]
				out = open(args.prefix+'.circular.'+str(l)+'.fasta','w')
				out.write(">%s_ciruclar_%s\n%s\n" %(s, str(l), sequence))
				out.close()
			
				#writing out file for testing circularity
				out = open(args.prefix+'.'+str(l)+'.for-testing.fasta','w')
				clip_from = int(clips[l][0][0])
				clip_to = int(clips[l][0][1])
				sequence = seqs[s][clip_to-500:clip_to]+seqs[s][clip_from:clip_from+500]
				out.write(">test_%s\n%s\n" %(str(l),sequence))
				out.close()
	else:
		print "\nDid not find any candidates in the specified length range\n"
