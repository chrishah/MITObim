#!/usr/bin/python

import sys,warnings
import argparse

VERSION="0.5"
DESCRIPTION='''

circules - checks for circularity in nucleotide sequences
version: v.%s

# Disclaimer: #
The script is currently in the beta phase. Any feedback is much appreciated.

''' %VERSION
kmers = []
dic = {}
seqs = {}
motifs_by_regions = {}
#auto=True


def write_clipped(seq_header, seq_seq, start, end, prefix):

	start = int(start)
	end = int(end)
	l = str(end - start)

	#writing out clipped sequence
	print "clip points: %i - %i - length: %s -> writing sequence to '%s.circular.%s.fasta'\n" %(start, end, l, prefix, l)
	sequence = seq_seq[start:end]
	out = open(prefix+'.circular.'+l+'.fasta','w')
	out.write(">%s_ciruclar_%s\n%s\n" %(seq_header, str(l), sequence))
	out.close()

	#writing out file for testing circularity
	out = open(prefix+'.'+l+'.for-testing.fasta','w')
	clip_from = start
	clip_to = end
	sequence = seq_seq[clip_to-500:clip_to]+seq_seq[clip_from:clip_from+500]
	out.write(">test_%s\n%s\n" %(str(l),sequence))
	out.close()

def roll_over(seq_seq, new_start, prefix):

	print "new start coordinate: %s -> writing sequence to '%s.roll.%s.fasta'\n" %(new_start, prefix, new_start)
	sequence = seq_seq[new_start-1:]+seq_seq[:new_start-1]
	out = open(prefix+'.rolled.'+str(new_start)+'.fasta', 'w')
	out.write(">%s_s%s_l%s\n%s\n" %(prefix, new_start, len(sequence), sequence))
        out.close()

parser = argparse.ArgumentParser(description=DESCRIPTION, prog='circules.py', 
	formatter_class=argparse.RawDescriptionHelpFormatter, 
	epilog='''   examples: 
	# check for circularity using a k-mer length of 31 - returns suggestion(s) for clipping points in putative circular sequences.
	circules.py -f test.fasta -k 31
	
	# check for circularity using k-mer lengths 31-51, stepsize 2bp - returns suggestion(s) for clipping points in putative circular sequences.
	circules.py -f test.fasta -k 31-51 -s 2
	
	# check for circularity using a k-mer length of 31. If length of suggested circular sequence is within +/- 1 percent of the expected length (15000 bp)
	# the clipped sequence will be written to a file called 'mytest.circular.fasta'. In addition a file 'mytest.for-testing.fasta' will be written. It contains
	# a 1000 bp sequence obtained by joining the first/last 500 bp of the proposed clipped circular sequence for additional evaluation if required.
	circules.py -f test.fasta -k 31 -l 15000 -p mytest

	# check for circularity using a k-mer length of 41. Extract if candidate is found in length range 9000 - 11000.
	circules.py -f test.fasta -k 41 -k 10000 --length_tolerance_percent 10 

	# clip sequence at specific clip points
	circules.py -f test.fasta -c 32,15430
	
	# roll circular sequence to new startposition, e.g. 46
	circules.py -f test.fasta -n 46

	''')

parser.add_argument("-f", "--fasta", help="fasta file containing the sequence to be evaluate.", metavar="<FILE>", action="store")
parser.add_argument("-k", "--kmer", help="kmer size. single number (default = 31) or range (e.g. 31-35).", type=str, default='31', action="store")
parser.add_argument("-s", "--kmer_step", help="kmer step size (default = 2).", type=int, default='2', action="store")
parser.add_argument("-p", "--prefix", help="prefix for output files (default = 'circular').", type=str, metavar="<STR>", action="store", default='circular')
parser.add_argument("-c", "--extract_by_coordinates", help="Coordinates for clipping of sequence in format 'startpos,endpos'. Clipped sequence will be written to file 'prefix.circular.fasta'.", type=str,  metavar="<INT,INT>", action="store")
parser.add_argument("-l", "--extract_by_length", help="expected length (in bp) of circular molecule. If a candidate of length expected (+-length tolerance if specified) is found, sequence will be clipped and written to file 'prefix.circular.fasta'.", type=int,  metavar="<INT>", default='0', action="store")
parser.add_argument("--length_tolerance_percent", help="percent length tolerance (e.g. 0.1, for 10 %%). Candidate fragments must have a length of 'expected length +/- t * expected length'. Default = 0.", type=float,  metavar="<FLOAT>", default='0', action="store")
parser.add_argument("--length_tolerance_absolute", help="absolute length tolerance (e.g. 1000). Candidate fragments must have a length of 'expected length +/- tolerance'. Default = 0.", type=int,  metavar="<INT>", default='0', action="store")
parser.add_argument("-n", "--newstart_roll", help="'roll' a putative ciruclar sequence to a specified new start point. Sequence will be written to file 'prefix.roll.{n}.fasta'.", type=int,  metavar="<INT>", default='0', action="store")
#parser.add_argument("-r", "--readpool", help="path to fastq reads to be used for evaluating circularity.", type=str, metavar="<FILE>", action="store")
#parser.add_argument("-v", "--verbose", help="turn verbose output on.", action="store_true")
parser.add_argument("--version", action="version", version=VERSION) 
args = parser.parse_args()

if len(sys.argv) < 2:	#if the script is called without any arguments display the usage
    print
    print "%s\n" %DESCRIPTION
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
		seqs[current]+=l.strip().upper()

if len(seqs) > 1:
	print "\nThe file '%s' contains %i sequences. The circules script is currently designed to only handle a single fasta sequence at a time - Sorry!\n" %(args.fasta, len(seqs))
	sys.exit()

print "%s\n" %DESCRIPTION

print "\n###########################\n\nProcessing sequence of length %s" %len(seqs[seqs.keys()[0]])

if args.newstart_roll:
	print "\n####################\n## Result output ##\n####################\n"
	print "\n'Rolling' circular sequence to new start coordinate .." 
	roll_over(seqs[seqs.keys()[0]], args.newstart_roll, args.prefix)
	print ""
	sys.exit()
	

if args.extract_by_coordinates:
	if not ',' in args.extract_by_coordinates:
		print "Expecting clipping coordinates delimited by comma (','), e.g. '-c 50,1567'"
		sys.exit()
	print "\n####################\n## Result output ##\n####################\n"
	print "\nExtracting sequence by coordinates .."
	clips = args.extract_by_coordinates.split(",")
	write_clipped(seqs.keys()[0], seqs[seqs.keys()[0]], clips[0], clips[1], args.prefix)
	print ""
	sys.exit()

for k in range(kmers[0],kmers[1],kmers[2]):
	clips = {}
	repeats = {}
	for s in seqs:
		print "\nCollecting %s-mers .." %k,
		i=0
		while i <= (len(seqs[s])-k):
			kmer = str(seqs[s][i:i+k])
#			print i,kmer,
			if not kmer[-1] in 'acgtACGT':
#				print "\tlast base of current kmer is ambiguous - jumping ahead by %s bases" %k
#				i+=k
				skip = k
				ok = False
				while ok == False:
					i+=skip
					kmer = str(seqs[s][i:i+k])
					skip = 0
#					print "%i %s" %(i, kmer)
					for pos in reversed(range(len(kmer))):
#						print pos,kmer[pos]
						if not kmer[pos] in 'acgtACGT':
#							print "%i\tambiguous" %pos
							skip = pos+1
							break
					if skip == 0:
						ok = True
				continue
			if not kmer in dic:
				dic[kmer] = []
			dic[kmer].append(i)
#			print "\tok"
			i+=1
		print "Done!"


	print "Evaluating %s-mers .." %k,

	for s in seqs:
		test = {}
		i=0
		while i <= (len(seqs[s])-k):
			kmer = str(seqs[s][i:i+k])
			
			if not kmer in dic:
				evaluate = False
			else:
				evaluate = True

			if evaluate:
				count = len(dic[kmer]) #this is the number of times the current kmer was found
				if count == 2:
					length = dic[kmer][1]-dic[kmer][0]
					if not length in clips:
						clips[length] = []
					clips[length].append(dic[kmer])

				elif count > 2:
					if not kmer in repeats:
						repeats[kmer] = []
					repeats[kmer].extend(dic[kmer])

			i+=1

		print "Done!\n"

	if repeats:

		print "\n###################\n## Repeat report ##\n###################"
		tandem = {}
		for kmer in repeats:
			repeats[kmer] = sorted(list(set(repeats[kmer])))
#			print "\n%s %s" %(kmer,repeats[kmer])

		#find blocks
			block = False
			blocks = [[],[]]
			dists = {}
			for i in range(len(repeats[kmer])-1):
				neighbour_dist = repeats[kmer][i+1] - repeats[kmer][i]
				if not neighbour_dist in dists:
					dists[neighbour_dist] = 0
				dists[neighbour_dist]+=1

#				print "%s vs %s" %(repeats[kmer][i],repeats[kmer][i+1])
				if neighbour_dist <= k:
					if not block:
						blocks[0].append(repeats[kmer][i])
						blocks[1].append(repeats[kmer][i+1]+k)
						block = True
					else:
						blocks[1][-1] = repeats[kmer][i+1]+k
				else:
					if not block:
						blocks[0].append(repeats[kmer][i])
						blocks[1].append(repeats[kmer][i]+k)
					
					block = False

			if not block:
				blocks[0].append(repeats[kmer][-1])
				blocks[1].append(repeats[kmer][-1]+k)


			tandem[int(blocks[0][0])] = blocks
			tandem[int(blocks[0][0])].append(dists)


		if tandem:
			motifs = []
			start_positions = sorted(tandem)

			active = False
			for i in range(len(start_positions)-1):
				
				if (start_positions[i+1] - start_positions[i]) == 1:
					if not active:
						motifs.append([start_positions[i],start_positions[i+1]])
						active = True
					else:
						motifs[-1][-1] = start_positions[i+1]

				else:
					if not active:
						motifs.append([start_positions[i],start_positions[i]])
					active=False
	
			if not active:
				motifs.append([start_positions[-1],start_positions[-1]])

			for m in motifs:
				
				block_start_end_coordinates = {}

				minstep = int(sorted(tandem[m[0]][2])[0])
				coordinates = range(m[0],m[1]+1)
				
				if minstep >= k:
					m[1]+=k-1
				elif minstep > (m[1]+1-m[0]):
					m[1]=m[0]+minstep-1

				for c in coordinates: #this loop finds the start and end position of each block
					for i in range(len(tandem[c][0])):
						if not i in block_start_end_coordinates:
							block_start_end_coordinates[i] = [tandem[c][0][i],tandem[c][1][i]]
						
						if tandem[c][0][i] < block_start_end_coordinates[i][0]:
							block_start_end_coordinates[i][0] = tandem[c][0][i]

						if tandem[c][1][i] > block_start_end_coordinates[i][1]:
							block_start_end_coordinates[i][1] = tandem[c][1][i]
						
				motif = seqs[s][m[0]:m[1]+1]
				if not motif in motifs_by_regions:
					motifs_by_regions[motif] = []
				print "\nmotif '%s', found multiple times in the following region(s) (start - end coordinate):" %seqs[s][m[0]:m[1]+1]
				

				for b in sorted(block_start_end_coordinates):
#					print b,block_start_end_coordinates[b]
					start = block_start_end_coordinates[b][0]
					end = block_start_end_coordinates[b][1]
					print "%i - %i (%i x)" %(start, end, (end-start)/(m[1]-m[0]+1))
#					print seqs[s][start:end+1]
#					print len(seqs[s])
					motifs_by_regions[motif].append([start,end])


	print "\n########################\n## Circularity report ##\n########################"
	nothing = True
	if len(clips) >=  1:
		print "\nFound %i candidate(s) for circularity supported by duplicated %i-mers:" %(len(clips), k)
		for l in sorted(clips):
			print "\t- suggested clip points %i - %i (length %s; supported by %i duplicted %i-mers); extract by adding '-l %i' or '-c %s,%s' to your command" %(int(clips[l][0][0]), int(clips[l][0][1]), l, len(clips[l])/2, k, l, clips[l][0][0], clips[l][0][1])
		nothing = False


	if motifs_by_regions:
		for m in motifs_by_regions:
			primes = [0]
			for i in range(len(motifs_by_regions[m])):
				if motifs_by_regions[m][i][0] == 0:
					primes[0] = i
				if motifs_by_regions[m][i][1] == len(seqs[s]):
					primes.append(i)
			if len(primes) == 2:
				print "\nFound motif '%s' at terminal positions (see 'Repeat report' above), which could indicate circularity." %m
				start_clip = motifs_by_regions[m][primes[0]][1]
				end_clip = motifs_by_regions[m][primes[1]][1]
				print "\t- suggested clip points: %s - %s (length %s); extract the clipped sequence by adding '-c %s,%s' to your command" %(start_clip, end_clip, end_clip-start_clip, start_clip, end_clip)
				print "\n\tThis will remove the repeated motif from the 5' end of the sequence.\n\tThe correct length can currently not be determined unambiguously because of the repeat. Longer reads could resolve this.\n"
				sys.exit()

	if nothing:
		print "\nDid not find any candidates for circularity using k = %s\n" %k
		sys.exit()
	
	print ""



if args.extract_by_length:
	print "\n####################\n## Result output ##\n####################\n"
	print "Extracting by expected length .."
	
	if args.length_tolerance_absolute:
		minlength = args.extract_by_length - args.length_tolerance_absolute
		maxlength = args.extract_by_length + args.length_tolerance_absolute
	elif args.length_tolerance_percent:
		minlength = float(args.extract_by_length)*(1-args.length_tolerance_percent)
		maxlength = float(args.extract_by_length)*(1+args.length_tolerance_percent)
	else:
		minlength = args.extract_by_length
		maxlength = minlength
	
	print "specified length: %s (tolerance: %s - %s) - Evaluating candidate lengths" %(args.extract_by_length, minlength, maxlength)
	for l in sorted(clips):
		if l <= maxlength and l >= minlength:
			pass
		else:
			del(clips[l])

	if clips:
		print "\nFound candidates in the expected length range"
		for s in seqs:
			for l in sorted(clips):
				write_clipped(s, seqs[s], clips[l][0][0], clips[l][0][1], args.prefix)

				#writing out clipped sequence
#				print "length: %s - clip points: %i - %i -> writing sequence to '%s.circular.%s.fasta'\n" %(str(l), int(clips[l][0][0]), int(clips[l][0][1]), args.prefix, str(l))
#				sequence = seqs[s][int(clips[l][0][0]):int(clips[l][0][1])]
#				out = open(args.prefix+'.circular.'+str(l)+'.fasta','w')
#				out.write(">%s_ciruclar_%s\n%s\n" %(s, str(l), sequence))
#				out.close()
			
#				#writing out file for testing circularity
#				out = open(args.prefix+'.'+str(l)+'.for-testing.fasta','w')
#				clip_from = int(clips[l][0][0])
#				clip_to = int(clips[l][0][1])
#				sequence = seqs[s][clip_to-500:clip_to]+seqs[s][clip_from:clip_from+500]
#				out.write(">test_%s\n%s\n" %(str(l),sequence))
#				out.close()
	else:
		print "\nDid not find any candidates in the specified length range\n"
