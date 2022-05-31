"""
Batch INDELible
Write control.txt file and execute indelible to simulate sequences.
"""
import sys
import os
from Bio import Phylo


# read tree from file
infile = sys.argv[1]

phy = Phylo.read(infile, 'newick')
for node in phy.get_nonterminals():
    node.name = None
tree_string = phy.format('newick')

scaling_factor = 1.  # average one substitution per site in tree

kappa1 = 4.0
kappa2 = 8.0

handle = open('control.txt', 'w')

# write minimal contents of INDELible control file
handle.write('[TYPE] NUCLEOTIDE 1\n')
handle.write('[SETTINGS]\n[output] FASTA\n[printrates] FALSE\n')
handle.write('[MODEL] model\n[submodel] TrN %f %f\n' % (kappa1, kappa2))

handle.write('[TREE] tree %s;\n' % tree_string.rstrip('0123456789.;:\n'))
handle.write('[treelength] %1.5f\n' % scaling_factor)
handle.write('[PARTITIONS] partitionname\n')
handle.write("  [tree model AY772699.txt]\n")
handle.write('[EVOLVE] partitionname 1 %s\n' % (infile, ))

handle.close()

#print 'running INDELIBLE'
os.system('indelible')
#os.system('/Users/art/src/INDELibleV1.03/src/indelible')  # can also use lldb
