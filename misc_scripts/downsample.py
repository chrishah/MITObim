#!/usr/bin/python


"""downsample

Author: Christoph Hahn (christoph.hahn@uni-graz.at)
February 2017

Extract a random subsample of ~ x % reads from fastq data.

The choice is based on a random number generator. For each fastq read, a random number between 1-100 will be generated. If the random number is smaller than the desired proportion in percent, the read will be kept, otherwise it will be discarded. So to extract ~15 % of the reads any read that gets a random number of <=15 will be kept, which will result in roughly 15% of the reads.

Subsamples can be taken from several fastq files at the same time. We allow to input paired end data in two separate files. If so specified subsamples will be taken so that the pairs will remain intact and the ouptut will be given in interleaved format.

Input fastq files can be compressed with gzipped. Mixed compressed / non-compressed input is possible except in the case of paired end data. In this case both read files need to be either compressed or non-compressed.

Examples:
# sample ~20 % of reads from three files
downsample.py -s 20 -r test.fastq.gz -r test2.fastq -r test3.fastq.gz > test.subsample_20.fastq
        
# sample ~30 % of reads from two files, and interleave reads from the two files on the fly
downsample.py -s 30 --interleave -r test_R1.fastq.gz -r test_R2.fastq.gz > test.interleaved.subsample_30.fastq

# sample ~40 % of reads from three files, defining a seed for the random number generator, to allow replication of the process.
downsample.py -s 20 --rand -421039 -r test.fastq.gz -r test2.fastq -r test3.fastq.gz > test.subsample_40.fastq

# sample ~20 % of reads from two files, compressing results on the fly.
downsample.py -s 20 -r test.fastq.gz -r test2.fastq | gzip > test.subsample_20.fastq.gz

"""

import sys
# import re
# import random

def parse_arguments():

	import sys
	import argparse

        VERSION="0.1"
        DESCRIPTION='''
        downsample.py - version: v.%s
        
        ''' %VERSION

        parser = argparse.ArgumentParser(description=DESCRIPTION, prog='downsample.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''   examples: 
        # sample ~20 % of reads from three files
        downsample.py -s 20 -r test.fastq.gz -r test2.fastq -r test3.fastq.gz
        
        # sample ~30 % of reads from two files, and interleave reads in output
        downsample.py -s 30 --interleave -r test_R1.fastq.gz -r test_R2.fastq.gz

	# sample ~40 % of reads from three files, defining a seed for the random number generator, to allow replication of the process.
	downsample.py -s 20 --rand -421039 -r test.fastq.gz -r test2.fastq -r test3.fastq.gz > test.subsample_40.fastq

	# sample ~20 % of reads from two files, compressing results on the fly.
	downsample.py -s 20 -r test.fastq.gz -r test2.fastq | gzip > test.subsample_20.fastq.gz

	#sample ~5% of reads from a single file that contains interleaved read data
	downsample.py -s 5 --interleave -r test.interleaved.fastq.gz | gzip > test.interleaved.subsample_5.fastq.gz

        ''')

        parser.add_argument("-r", "--reads", help="Readfile (flag can be used repeatadly to process several files", metavar="<FILE>", action="append")
        parser.add_argument("-s", "--sample", help="Desired size of subsample in percent (1-100; default = 50)", type=int, metavar="<INT>", action="store", default=50)
        parser.add_argument("--interleave", help="Optional. In case of two input files, data will be interleaved from these in the output. Otherwise data will be treated as already interleaved.", action="store_true")
        parser.add_argument("--seed", help="Optional. Seed for random number generator", metavar="<INT>", type=int,  action="store")

        parser.add_argument("--version", action="version", version=VERSION)

	if not parser.parse_args().reads or len(sys.argv) == 1:
                print
                parser.print_usage()
                print
                sys.exit(1)

	return parser.parse_args()

def check_args(args):

        if args.sample < 1 or args.sample > 100:
              sys.exit("\n only sample size 1-100 is valid\n")


def set_seed(seed):

	import random

	if not seed:
		seed = random.randint(-100000000,100000000)

	sys.stderr.write("seed for random number generator is: %i\n" %seed)
        random.seed(seed)
			
def decide(string, percent):
        
	import random

	if (random.randint(1,100) <= percent):
                print string,

def sample_interleave(file1, file2, percent):

    import re

    while True:
        out = ''
        line = file1.readline()
        if line.strip() == "":
            break
        out+=re.sub(r" 1:n.*", "/1",line)

        for i in xrange(3):
            out+=re.sub(r" 2:n.*","/2",file1.readline())

        for i in xrange(4):
            out+=re.sub(r" 2:n.*","/2",file2.readline())

        decide(out, percent)

def sample(fi, percent, step):

    import re
    
    while True:
        out = ''
        line = fi.readline()
        if line.strip() == "":
            break
        out+=line
        for i in xrange(step-1):
                out+=fi.readline()

        decide(out, percent)

def main():

	import sys

	args = parse_arguments()
	check_args(args)

        sys.stderr.write("\ndownsampling to %i percent\n" %args.sample)

	set_seed(args.seed)

	if args.interleave and len(args.reads) == 2:
                sys.stderr.write("interleaving sample from input files %s and %s\n" %(args.reads[0], args.reads[1]))
                if args.reads[0][-2:] == "gz":
                        import gzip
                        with gzip.open(args.reads[0]) as f1:
                                with gzip.open(args.reads[1]) as f2:
                                        sample_interleave(f1, f2, args.sample)
                else:
                        with open(args.reads[0]) as f1:
                                with open(args.reads[1]) as f2:
                                        sample_interleave(f1, f2, args.sample)

                f1.close()
                f2.close()

	else: #that is all other cases
    		if args.interleave:
			sys.stderr.write("you indicated interleaved input file(s) -> stepsize = 8 lines\n")
			step = 8
    		else:
			sys.stderr.write("you indicated single end data -> stepsize = 4 lines\n")
			step = 4

                for readsfile in args.reads:
                        if readsfile[-2:] == "gz":
                                import gzip
                                f = gzip.open(readsfile)
                        else:
                                f = open(readsfile)

                        sample(f, args.sample, step)
                        f.close()

	sys.stderr.write("Done!\n\n")

if __name__ == '__main__':
        sys.exit(main())

