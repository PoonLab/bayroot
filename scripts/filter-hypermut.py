from glob import glob
from poplars import hypermut
from poplars.common import *

# multiple alignments (MAFFT)
# for f in `ls ZM*.fa`; do mafft --reorder $f > $f.maf; done
for f in glob('data/*.maf'):
    print(f)
    with open(f) as handle:
        fasta = convert_fasta(handle)

    # generate consensus sequence from seroconversion sequences
    sc = [[h, s] for h, s in fasta if '_SC_' in h]
    conseq = consensus(sc)

    # filter data for hypermutated sequences
    results = hypermut.hypermut(f, cons=conseq)
    results = dict([(row['seq_name'], row) for row in results])
    print(len(results))
    filtered = [(h, s) for h, s in fasta if results[h]['p_value'] > 0.05]
    print(len(filtered))

    outfile = open(f.replace('.fa.maf', '.fa.hyp'), 'w')
    for h, s in filtered:
        outfile.write('>{}\n{}\n'.format(h, s))
    outfile.close()

# for f in `ls data/ZM*.fa.hyp`; do iqtree -m GTR+G -s $f; done
