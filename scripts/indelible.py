import sys
import os
from Bio import Phylo
import tempfile
import argparse
import subprocess

# command line interface
parser = argparse.ArgumentParser("Simulate sequences using INDELible with a "
                                 "user-specified tree and root sequence.")
parser.add_argument("infile", type=str, help="Path to file containing Newick "
                    "tree string")
parser.add_argument("--root", type=str, default="AY772699.txt",
                    help="Path to plain text file containing sequence "
                         "to assign to root.")
parser.add_argument("--sf", type=float, default=1.0,
                    help="Scaling factor: expected number "
                    "of substitutions per site across entire tree length. "
                    "Defaults to 1.")
parser.add_argument("--k1", type=float, default=4.0,
                    help="Rate bias for C-T transitions.")
parser.add_argument("--k2", type=float, default=8.0,
                    help="Rate bias for A-G transitions.")
args = parser.parse_args()

# read tree from file
phy = Phylo.read(args.infile, 'newick')
for node in phy.get_nonterminals():
    node.name = None
tree_string = phy.format('newick').rstrip('0123456789.;:\n')

template = f"""[TYPE] NUCLEOTIDE 1
[SETTINGS]
  [output] FASTA
  [printrates] FALSE
[MODEL] model
  [submodel] TrN {args.k1} {args.k2}
  [statefreq] 0.2 0.2 0.4 0.2
[TREE] tree {tree_string};
  [treelength] {args.sf}
[PARTITIONS] pname
  [tree model {args.root}]
[EVOLVE] pname 1 {args.infile}
"""

# write control text to a temporary file
tmpfile = tempfile.NamedTemporaryFile('wt', delete=False)
tmpfile.write(template)
tmpfile.close()

# run INDELible with this control file
subprocess.check_call(['indelible', tmpfile.name])
