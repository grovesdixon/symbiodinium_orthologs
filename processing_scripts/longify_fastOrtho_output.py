#!/usr/bin/env python

#import modules
import argparse





##############################
###### DEFINE FUNCTIONS ######
##############################

def read_fastOrtho(orthoFile):
    '''Read through the fastortho .end output and output the sequences for each orthogroup
    '''
    print('\nParsing fastOrtho file {}...'.format(orthoFile))
    totalGroups = 0
    totalSeqs = 0

    #open the fastOrtho groups file
    with open(orthoFile, 'r') as infile:
        with open(outfile,  'w') as out:
            out.write("orthoGroup\tcontig")
            #each line contains one orthologous group
            for line in infile:
                orthoDat = line.strip("\n").split(":")[1]
                oList = orthoDat.split() #the list of the sequences in this orthologous group (with appended speceis names)
                line = line.strip("\n").split()
                groupName = line[0]
                numGenes = int(line[1].strip("("))
                numTaxa = int(line[2].strip("genes,"))
                proteinSeqs = []
                totalGroups += 1
                totalSeqs += numGenes
                for homolog in oList:
                    seqid = homolog.split("(")[0]
                    out.write('\n{}\t{}'.format(groupName, seqid))
    print('Total orthologous groups = {}'.format(totalGroups))
    print('Total sequences = {}'.format(totalSeqs))
    print('Results saved to {}\n'.format(outfile))




    
##################################
############## MAIN ##############
##################################

if __name__ == '__main__':

    #START RUN TIME CLOCK

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Check a tree to ensure that a given set of taxa are monophyletic
    '''
    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, dest = 'input', help = 'ortholog file from FastOrtho usually ends with .end')
    parser.add_argument('-o', required = True, dest = 'output', help = 'name for output file')
    
    #----- PARSE ARGUMENTS -----#
    args = parser.parse_args()
    infile = args.input
    outfile = args.output

    #------ RUN FUNCTIONS ------#
    read_fastOrtho(infile)
            


