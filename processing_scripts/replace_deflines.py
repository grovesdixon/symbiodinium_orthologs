#!/usr/bin/env python
##replace_deflines.py
##last updated 12-15-19
    ###revised to write out as iterator rather than compiling list
##Groves Dixon

#import modules
import time
import argparse
from sys import argv
from sys import exit
from Bio import SeqIO



##############################
###### DEFINE FUNCTIONS ######
##############################



##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Replace the definition lines in a fasta file with a prefix and increasing numbers.
    Optionally remove description with -clean True
    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-fa', required = True, dest = 'input_file', help = 'The the input fasta')
    parser.add_argument('-prefix', required = True, dest = 'prefix', help = 'The desired prefix.')
    parser.add_argument('-o', required = True, dest = 'output_file', help = 'The desired output name.')
    parser.add_argument('-clean', required = False, default=False, dest = 'cleanLines', help = 'Change to "yes" to remove all other information from sequence description')

    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    infileName = args.input_file
    prefix = args.prefix
    outfileName = args.output_file
    cleanLines = args.cleanLines
    matchTable = "{}_{}_seqTable.tsv".format(infileName, outfileName)


    #--- RUN ---#
    pairList = []
    seqNum = 0
    fasSeqs = SeqIO.parse(open(infileName), 'fasta')
    with open(outfileName, 'w') as out:
        for seq in fasSeqs:
            seqNum += 1
            seqId = seq.id
            newId = "{}_{}".format(prefix, seqNum)
            pairList.append("{}\t{}".format(newId, seqId))
            seq.id = newId
            if cleanLines:
                seq.description = ''
            outRec = seq.format("fasta-2line")
            out.write(outRec)

    with open(matchTable, 'w') as out:
        for pair in pairList:
            out.write(pair + "\n")



