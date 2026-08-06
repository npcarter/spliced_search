#! python

import argparse
import os

parser = argparse.ArgumentParser(description='Extracts the primary and single-sequence sequences from a file of training data')
parser.add_argument('infile')
parser.add_argument('outfile')
args = parser.parse_args()

counter = 0
writeThis = False

with open (args.infile, 'r') as infile, open(args.outfile, 'a') as outfile:
    for line in infile:
        if counter % 3 == 0:
            fields = line.split('^')
            if fields[-1].strip() == 'primary' or fields[-1].strip() == 'single':
                writeThis = True
            else:
                writeThis = False
        if writeThis:
            outfile.write(line)
        counter +=1 