#!/usr/bin/python
# encoding:utf8
# authors: Erik Garrison, Sébastien Boisvert
# modified by github@cypridina on 20151104 to work with MITObim
"""This script takes two fastq or fastq.gz files and interleaves them
Usage:
    interleave-fasta fasta_file1 fasta_file2
"""

import sys,re 

def interleave(f1, f2):
    """Interleaves two (open) fastq files.
    """
    while True:
        line = f1.readline()
        if line.strip() == "":
            break
        print re.sub(r" 1:N.*", "/1",line.strip())
        
        for i in xrange(3):
            print re.sub(r" 2:N.*","/2",f1.readline().strip())
        
        for i in xrange(4):
            print re.sub(r" 2:N.*","/2",f2.readline().strip())

if __name__ == '__main__':
    try:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
    except:
        print __doc__
        sys.exit(1)

    if file1[-2:] == "gz":
        import gzip
        with gzip.open(file1) as f1:
            with gzip.open(file2) as f2:
                interleave(f1, f2)
    else:
        with open(file1) as f1:
            with open(file2) as f2:
                interleave(f1, f2)
    f1.close()
    f2.close() 
