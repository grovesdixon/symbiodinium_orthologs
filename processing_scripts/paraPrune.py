#!/usr/bin/env python3
##paraPrune.py
##last upadated 3-9-20
    ##updated for python3, and to use pandas to output a dataframe of 
     #collapsed contigs for summing later in RNAseq analysis
    ##12-9-19
        ### added two more steps of pruning with prune_any_duplicate() and pull_perfect_nodes()
    ##9-25-18
        ###added function save_tree_subset() for spot checking results
        ##removed some junk from the script
        ##Also caught that I was not identifying duplcated species within polytomies and fixed
##Groves Dixon


#import modules
import time
import argparse
from sys import argv
from sys import exit
from Bio import Phylo
import numpy as np
import os
import datetime
import pandas as pd





##############################
###### DEFINE FUNCTIONS ######
##############################

def read_lengths():
    """Read in the length file. This is just for choosing longest paralogs"""
    ldict = {}
    with open(lengthsInput, 'r') as infile:
        for line in infile:
            line=line.strip("\n").split()
            ldict[line[0]] = float(line[1])
    return ldict

def unique(list1):
    x = np.array(list1)
    return(list(np.unique(x)))

def read_in_trees(treeFileList):
    print("\nReading in trees...")
    emptyFileList = []
    treeList = []
    for treeFile in treeFileList:
        try:
            tree = Phylo.read(treeFile, 'newick')
        except ValueError:
            emptyFileList.append(treeFile)
            print("Warning, tree file {} doesn't seem to have a tree.".format(treeFile))
            continue
        treeList.append(tree)
        tree.name = treeFile.split(".newick")[0]
    print("\t{} tree files were empty, suggest checking these:".format(len(emptyFileList)))
    for t in emptyFileList:
        print(t)
    print("\t{} total trees for analysis".format(len(treeList)))
    return treeList


def save_tree_subset(treeList, treeType, subDirName):
    """Function to output a subset of trees for reference"""
    print("\tWriting out {} {} trees to subdirectory {}...".format(len(treeList), treeType, subDirName))
    if not os.path.exists(subDirName):
        os.makedirs(subDirName)
    with open("{}/{}.ascii".format(subDirName, subDirName), 'w') as asciiOut:
        for treeToWrite in treeList:
            Phylo.write(treeToWrite, "{}/{}_paraPruned.newick".format(subDirName, treeToWrite.name), 'newick')
            asciiOut.write("============\n")
            asciiOut.write("{}:\n".format(treeToWrite.name))
            Phylo.draw_ascii(treeToWrite, file=asciiOut)



def find_single_copy(treeList):
    singleCopyList = []
    hasParaList = []
    for tree in treeList:
        terminalSeqList = [term.name for term in tree.get_terminals()]
        terminalSppList = [term.name.split("_")[0] for term in tree.get_terminals()]
        uSppList = unique(terminalSppList)
        singleCopy = len(uSppList) == len(terminalSppList)
        if singleCopy:
            singleCopyList.append(tree)
        else:
            hasParaList.append(tree)
        # print("-----------")
        # print(tree.name)
        # print(terminalSeqList)
        # print(terminalSppList)
        # print(uSppList)
        # print("Single Copy = {}".format(singleCopy))


    print("\t{} total trees checked for paralogs".format(len(treeList)))
    print("\t{} trees were single copy".format(len(singleCopyList)))
    print("\t{} got paralogs".format(len(hasParaList)))
    return(singleCopyList, hasParaList)


def get_longest_seq(leafList):
    """use sequence lengths from read_lengths() to pick longest
    version for an single-species clade"""
    nameList = [leaf.name for leaf in leafList]     #get the leaf names, which are the sequence names
    lengthList = [ldict[name] for name in nameList] #get the sequence lengths pulled with read_lengths()
    for i in reversed(range(len(lengthList))):      #reversed so first becomes assigned if they are equal
        if lengthList[i] == max(lengthList):
            longestSeq = leafList[i]
    # print("------------")
    # print(leafList)
    # print(nameList)
    # print(lengthList)
    # print("longest = {}".format(longestSeq.name))
    return(longestSeq)


def find_dup_in_polytomy(nodeTerminals, nodeSpp, uNodeTaxa):
    """Function to identify leaves to remove when species are duplicated in polytomires.
    Returns a set of leafs to remove, along with the longest contig from that species
    that was kept. These can then be added to the collapsed dataframe in collapse_monophyletic()"""
    toRemoveList = []
    keepList = []
    tdict = {}
    for unt in uNodeTaxa:
        tdict[unt] = []
        for ndx in range(len(nodeSpp)):
            if nodeSpp[ndx]==unt:
                tdict[unt].append(nodeTerminals[ndx])
    for unt in uNodeTaxa:
        termList = tdict[unt]
        if len(termList) > 1: #ie if there's more than one contig for this species
            toKeep = get_longest_seq(termList)
            termList.remove(toKeep)
            toRemoveList.append(termList)
            keepList.append(toKeep)
    return(toRemoveList, keepList)





def collapse_monophyletic(treeList):
    print("\nLooking for single-species clades...")
    collapseCount = 0
    collapsedList = []
    collapsedDf = pd.DataFrame(columns = ['orthogroup', 'kept', 'collapsed']) #dataframe to store collapsed contigs for getting back RNAseq counts

    for tree in treeList:
        nodeList = tree.get_nonterminals()
        for node in nodeList:
            nodeTerminals = node.get_terminals()
            nodeSpp = [leaf.name.split("_")[0] for leaf in nodeTerminals]
            uNodeTaxa = unique(nodeSpp)

            # print('--')
            # print("N terminals = {}".format(len(nodeTerminals)))
            # print("N unique = {}".format(len(uNodeTaxa)))
            # print("terminal = {}".format(node.is_preterminal()))


            #use the unique set of taxa in this node as a test for a single-species clade
            if len(uNodeTaxa) == 1:

                # print("===============================")
                # print("Found one!")
                # Phylo.draw_ascii(tree)

                nameList = [leaf.name for leaf in nodeTerminals]
                toKeep = get_longest_seq(nodeTerminals)
                toRemove = nodeTerminals
                toRemove.remove(toKeep)

                for leaf in toRemove:
                    tree.prune(leaf)
                    collapsed_data = pd.DataFrame({
                    'orthogroup':[tree.name],
                    'kept':[toKeep.name],
                    'collapsed':[leaf.name]})
                    # print('collapsed data:')
                    # print(collapsed_data)
                    collapsedDf = collapsedDf.append(collapsed_data)
                collapseCount += 1
                # print('------------------------')
                # Phylo.draw_ascii(tree)

            #solve possible polytomy problem:
            elif node.is_preterminal() and len(nodeSpp) > len(uNodeTaxa):         
                toRemove, toKeep = find_dup_in_polytomy(nodeTerminals, nodeSpp, uNodeTaxa)

                # print("==================")
                # print("Found a polytomy with duplicated species!")
                # Phylo.draw_ascii(tree)
                # print(nodeTerminals)
                # print(nodeSpp)
                # print("These leafs need to be removed:")
                # print(toRemove)
                # print("These leafs to be kept:")
                # print(toKeep)
                
                #do the pruning
                for trml in toRemove:
                    for leaf in trml:
                        tree.prune(leaf)
                        collapseCount += 1

                #record contigs collapsed from polytomy
                for i in range(len(toRemove)):
                    collapsed_spp_nodes = toRemove[i] #sublist of contigs for this species (probably just one node)
                    kept_node = toKeep[i]
                    for collapsed_node in collapsed_spp_nodes:
                        collapsed_data = pd.DataFrame({
                        'orthogroup':[tree.name],
                        'kept':[kept_node.name],
                        'collapsed':[collapsed_node.name]})
                        # print('collapsed data:')
                        # print(collapsed_data)
                        collapsedDf = collapsedDf.append(collapsed_data)



            else:
                pass

        collapsedList.append(tree)
    print("\t{} single-species clades were collapsed into single sequence".format(collapseCount))
    return(collapsedList, collapsedDf)

def output_collapsed(collapsed_df):
    """Function to output a dataframe of collapsed contigs.
    Idea here is that counts for these in an RNAseq counts matrix
    should be counted for the expression levels, since they are being
    contidered the 'same' gene"""
    print('\nOutputting collapsed contigs. These can be summed back in later in RNAseq analysis.')
    print('Total collapsed contigs = {}'.format(collapsed_df.shape[0]))
    print('Writing these out to collapsed_contigs.tsv...')
    collapsed_df.to_csv("collapsed_contigs.tsv", sep="\t", index=False)


def prune_uninportant(treeList):
    print("\nPruning away anemone species and Adig references...")
    unimportant = ["Aelegantis", "Apallida", "Nvectensis", "XP"]
    prunedUnimport = 0
    prunedPolytomy = 0
    revisedTrees = []
    for tree in treeList:
        originalTree=tree
        leafList = tree.get_terminals()
        for leaf in leafList:
            if leaf.name.split("_")[0] in unimportant:
                tree.prune(leaf)
                prunedUnimport += 1
        revisedTrees.append(tree)

    #removal of these branches can result in polytomies
    #these can then have multiple instances of single species in them
    #so remove those as if they were single-species clades as in collapse_monophyletic()
    for tree in treeList:
        nodeList = tree.get_nonterminals()
        for node in nodeList:
            if node.is_preterminal():
                leafList = node.get_terminals()
                leafSppList = [leaf.name.split("_")[0] for leaf in leafList]
                uSppList = unique(leafSppList)

                #check if any species appears more than once in the polytomy
                if len(leafSppList) > len(uSppList):
                    # print("FOUND POLYTOMOY WITH REPEATED SPECIES")
                    # Phylo.draw_ascii(tree)

                    sppCountList = [leafSppList.count(spp) for spp in leafSppList]
                    repeatSppLeafs = []
                    for i in range(len(sppCountList)):
                        if sppCountList[i] > 1:
                            repeatSppLeafs.append(leafList[i])

                    toKeep = get_longest_seq(repeatSppLeafs)
                    toRemove = repeatSppLeafs
                    toRemove.remove(toKeep)
                    # print("LEAFS TO REMOVE:")
                    # print(uSppList)
                    # print(leafSppList)
                    # print(sppCountList)
                    # print(repeatSppLeafs)
                    # print(toRemove)
                    for leaf in toRemove:
                        tree.prune(leaf)
                        prunedPolytomy += 1
                    # print("---------------")
                    # Phylo.draw_ascii(tree)
                else:
                    continue
            else:
                continue
    print("\t{} leafs from unimportant group pruned away".format(prunedUnimport))
    print("\t{} repeated species within resulting polytomies were pruned".format(prunedPolytomy))
    return(revisedTrees)



def remove_poorly_supported(treeList):
    print("\nPruning away repeated species branches that are not well supported...")
    for tree in treeList:
        nodeList = tree.get_nonterminals()
        print(nodeList)
        for node in nodeList:
            if node.is_preterminal():
                leafList = node.get_terminals()
                leafSppList = [leaf.name.split("_")[0] for leaf in leafList]
                uSppList = unique(leafSppList)
                print(node)


def prune_any_duplicate(treeList):
    """Function to prune away ALL species with repeated sequences if repeated seqeunces are relatively rare in tree.
    Here relatively rare means that the median count of sequences rounded to nearest integer is less than 2.
    Otherwise tree is assumed to contain multiple genes, which can maybe be pulled out with pull_clades.py"""
    print("\nPruning away ANY repeated species branches in trees with median number of appearances == 1...")
    revisedTrees = []
    prunedAnyDups = 0
    manySeqTrees = 0
    for tree in treeList:
        originalTree=tree
        leafList = tree.get_terminals()
        seqNameList = [leaf.name for leaf in leafList]
        sppList = [x.split('_')[0] for x in seqNameList]
        sppArr = np.array(sppList)
        allSpp, counts = np.unique(sppArr, return_counts=True)
        allSpp.sort()
        totSpp = len(allSpp)
        med = int(np.median(counts))
        if med < 2:
            prunedAnyDups += 1
            dupSpecies = []
            for i in range(len(allSpp)):
                c=counts[i]
                if c > 1:
                    dupSpecies.append(allSpp[i])
            for l in leafList:
                sppName = l.name.split('_')[0]
                if sppName in dupSpecies:
                    tree.prune(l)
        else:
            pass
        revisedTrees.append(tree)
    return(revisedTrees)


def parse_node(node):
  """get the species, unique species, and count for a node"""
  nodeTerminals = node.get_terminals()
  nodeSpp = np.array([leaf.name.split("_")[0] for leaf in nodeTerminals])
  uSpp, uCounts = np.unique(nodeSpp, return_counts=True)
  nSpp = len(uSpp)
  nodeSppList = list(nodeSpp)
  nodeSppList.sort()
  branchLengths = [t.branch_length for t in nodeTerminals]
  leafNames = [t.name for t in nodeTerminals]
  return(nodeSppList, uSpp, nSpp, branchLengths, leafNames)


def read_in_perfect_node_trees(perfectNodeList, subDirName):
    """Function to read back in subTrees identified in pull_perfect_nodes()"""
    perfectNodeTrees = []
    for treeFile in perfectNodeList:
        tree = Phylo.read(treeFile, 'newick')
        newName = treeFile.split('{}/'.format(subDirName))[1].split('_paraPruned.newick')[0]
        tree.name = newName
        perfectNodeTrees.append(tree)
    return(perfectNodeTrees)



def pull_perfect_nodes(treeList, subDirName):
    """Function to pull any perfect nodes from the tree,
    where perfect means the node has exactly 1 member from all species in the tree"""
    print("\nPulling out any perfect subtrees from remaining trees with paralogs")
    print('\tChecking through {} trees with paralogs'.format(len(treeList)))
    if not os.path.exists(subDirName):
        os.makedirs(subDirName)
    with open("{}/{}.ascii".format(subDirName, subDirName), 'w') as asciiOut:
        remainingTreeList = []
        perfectNodeList = []
        perfectNodeTrees = 0
        totalPerfectNodes = 0
        completelyPerfect = 0
        for tree in treeList:
            treeName = tree.name
            treePerfectNodes = 0
            # print('==================')
            # print(treeName)
            # Phylo.draw_ascii(tree)
            hadPerfect = False
            hasLeavesLeft = True
            leafList = tree.get_terminals()
            seqNameList = [leaf.name for leaf in leafList]
            sppList = [x.split('_')[0] for x in seqNameList]
            sppArr = np.array(sppList)
            allSpp, counts = np.unique(sppArr, return_counts=True)
            allSpp.sort()
            totSpp = len(allSpp)
            nodeList = tree.get_nonterminals()
            subTreeList = []
            for node in nodeList:
                nodeSppList, uSpp, nSpp, branchLengths, leafNames = parse_node(node)
                if nodeSppList == list(allSpp):
                    hadPerfect = True
                    treePerfectNodes += 1
                    totalPerfectNodes += 1
                    subTreeOutName = "{}/{}.{}_paraPruned.newick".format(subDirName, treeName, treePerfectNodes)
                    # print('-----')
                    # print(node.get_terminals())
                    # print(nodeSppList)
                    # print(node)
                    Phylo.write(node, subTreeOutName, 'newick')
                    asciiOut.write("============\n")
                    asciiOut.write("{}:\n".format(subTreeOutName))
                    if sum(branchLengths) > 0.0:    #had to add this in because you can't draw a zero length tree
                        Phylo.draw_ascii(node, file=asciiOut)
                    else:
                        asciiOut.write("-".join(leafNames) + "\n")
                    perfectNodeList.append(subTreeOutName)
                    leafsLeft = tree.get_terminals()
                    #remove the perfect subtree by pruning its leaves
                    if len(leafsLeft)==len(leafNames): #added this for case when tree is fully pruned into perfect clades
                        hasLeavesLeft = False
                        completelyPerfect += 1
                        continue
                    else:
                        for l in node.get_terminals():
                            tree.prune(l)
            perfectNodeTrees += int(hadPerfect)
            if hasLeavesLeft:
                leafList = tree.get_terminals()
                if len(leafList) > 3:
                    remainingTreeList.append(tree)

        print('\t{} trees had at least one perfect subtree that was made into its own ortholog'.format(perfectNodeTrees))
        print('\t{} total perfect nodes were found'.format(totalPerfectNodes))
        print('\t{} total trees were separated completely into perfect subtrees'.format(completelyPerfect))
        print('\t{} total perfect subtrees found among {} full trees'.format(len(perfectNodeList), perfectNodeTrees))
        print('\t{} trees still had at least 3 terminal branches left after pulling perfect subtrees'.format(len(remainingTreeList)))
        print('\t-------------')
        print('\tfinal check for single copy trees')

    #read back in the perfect nodes as trees
    perfectNodeTrees = read_in_perfect_node_trees(perfectNodeList, subDirName)
    return(perfectNodeTrees, remainingTreeList)


def output_revised_groups(treeList):
    outName = 'singleCopyOrthos.txt'
    print("\nSingle copy orthos saved as {}...".format(outName))
    with open(outName, 'w') as out:
        for tree in treeList:
            leafList = tree.get_terminals()
            seqNameList = [leaf.name for leaf in leafList]
            for sn in seqNameList:
                out.write("{}\t{}\n".format(tree.name, sn))
        



##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
    Start_time = time.time() 

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Look through gene trees made from protein alignments of sequences in orthologous groups.
    Output trees and sequence sets for single-copy orthologs.
    Then collapse clades that are monophyletic for a given species and keep only the longest sequence to get more.
    At each step, save the trees in subdirectories.
    
    Outputs:
    singleCopyOrthos.txt - conatains the final set of single copy orthos after pruning


    '''
    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-trees', required = True, dest = 'tree_files', nargs="+", help = 'Glob to the orthogroup trees')
    parser.add_argument('-l', required = True, dest = 'length_file', help = 'A table of lengths for all the OFRs used in the blast for FastOrtho')
    parser.add_argument('-subsets', required = False, default=False, dest = 'output_subset_trees', help = 'Boolean for whether to make copies of trees from each category during pruning process. Does not effect final result, just whether there are lots of copied files for reference.')
    
    #----- PARSE ARGUMENTS -----#
    args = parser.parse_args()
    treeFileList = args.tree_files
    lengthsInput = args.length_file
    makeSubsets = args.output_subset_trees
    # overwrite = args.overwrite


    #------ RUN FUNCTIONS ------#

    #initialize
    print("\nRunning paraPrune.py...")
    if makeSubsets:
        print("\nWill save subsets of trees at each step.")

    #read in trees
    treeList0 = read_in_trees(treeFileList)


    #read in trees and get intial single-copy orthogroups
    print("\nGathering single copy trees...")
    rawSingleCopy, remainingTreeList = find_single_copy(treeList0)

    #save the trees that were initially single copy or had paralogs
    if makeSubsets:
        save_tree_subset(rawSingleCopy, "initally single copy", "A_init_single_copy_trees")
        save_tree_subset(remainingTreeList, "initally with paralogs", "A_init_paralog_trees")

    
    #collapse single-species clades, keeping only longest sequence
    #(assume that these are isoforms that failed to cluster in cd-hit)
    ldict = read_lengths()
    collapsedTreeList, collapsed_df = collapse_monophyletic(remainingTreeList) #does the collapsing
    output_collapsed(collapsed_df)

    #now recheck the trees for single copies after the collapsing
    print("\nGathering additional single copy trees after collapsing single-species nodes...")
    collapseSingleCopy, remainingTreeList = find_single_copy(collapsedTreeList)
    allSingleCopy = rawSingleCopy + collapseSingleCopy
    print("\t{} total single copy orthogroups".format(len(allSingleCopy)))

    #save trees that were single copy after pruning, and still had paralogs
    if makeSubsets:
        save_tree_subset(collapseSingleCopy, "single-copy post pruning", "B_post_prune_single_copy_trees")
        save_tree_subset(remainingTreeList, "with paralogs post pruning", "B_post_prune_paralog_trees")


    #prune unimportant branches
    prunedTrees = prune_uninportant(remainingTreeList)
    print("\nGathering additional single copy trees after pruning anemones and Adig reference...")
    postPruneSingleCopy, remainingTreeList = find_single_copy(prunedTrees)
    allSingleCopy = rawSingleCopy + collapseSingleCopy + postPruneSingleCopy
    print("\t{} total single copy orthogroups".format(len(allSingleCopy)))

    #save trees that were single copy after pruning, and still had paralogs
    if makeSubsets:
        save_tree_subset(postPruneSingleCopy, "single-copy post pruning unimportant", "C_post_prune_single_copy_trees")
        save_tree_subset(remainingTreeList, "with paralogs post pruning", "C_post_prune_paralog_trees")

    #prune away small stray multi-sequence species
    allDupPrunedTrees = prune_any_duplicate(remainingTreeList)
    print("\nGathering additional single copy trees after pruning ANY remainging duplicates if they are rare...")
    postAnySingleCopy, remainingTreeList = find_single_copy(allDupPrunedTrees)
    allSingleCopy = rawSingleCopy + collapseSingleCopy + postPruneSingleCopy + postAnySingleCopy
    print("\t{} total single copy orthogroups".format(len(allSingleCopy)))

    #save trees that were single copy after pruning, and still had paralogs
    if makeSubsets:
        save_tree_subset(postAnySingleCopy, "single-copy post pruning ANY", "D_post_prune_single_copy_trees")
        save_tree_subset(remainingTreeList, "with paralogs post pruning ANY", "D_post_prune_paralog_trees")

    #pull out any 'perfect' subtrees, which have 1 sequence from all species in the tree
    perfectNodeTrees, leftOver = pull_perfect_nodes(remainingTreeList, "E_perfect_node_subtrees")
    postPerfectSingleCopy, remainingTreeList = find_single_copy(leftOver) #see if any single copies arose after pull perfect nodes
    allSingleCopy = rawSingleCopy + collapseSingleCopy + postPruneSingleCopy + postAnySingleCopy + perfectNodeTrees + postPerfectSingleCopy
    print("\t{} total single copy orthogroups".format(len(allSingleCopy)))

    #save trees that were single copy after pruning, and still had paralogs
    if makeSubsets:
        save_tree_subset(postPerfectSingleCopy, "single-copy post pull perfect", "E_post_perfect_node_single_copy_trees")
        save_tree_subset(remainingTreeList, "with paralogs post pruning ANY", "F_paralog_trees")


    #output final results
    output_revised_groups(allSingleCopy)

    #return time to run
    print("\nDone.")
    Time = time.time() - Start_time
    print('\nTime took to run: {}\n'.format(Time))




