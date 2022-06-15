from datetime import datetime, date
import sys

# ZMN133M_16Nov2011_2A_9_E
# ZM = country code (Zambia)
# N133M = patient code (M/F=gender)
# 16Nov2011 = sample collection date
# 2A - ART2? PBMC, post-ART
# SC = seroconversion sample (plasma),
# TF = inferred transmitted founder, consensus of sequences from earliest sample
# CP = plasma, drug naive (baseline?)
# CC = PBMC, drug naive
# 1Y = year 1, plasma, drug naive

def convert_fasta(handle):
    result = []
    sequence = ''
    for line in handle:
        if line.startswith('$'):  # skip header line
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h, sequence])
                sequence = ''  # reset
            h = line.strip('>#\n')
        else:
            sequence += line.strip('\n').upper()

    result.append([h, sequence])  # handle last entry
    return result

with open('data/zambia.fa') as handle:
    fasta = convert_fasta(handle)

# separate by patient
by_patient = {}
for h, s in fasta:
    accno, _, _, qname = h.split()[:4]
    if qname.endswith('_TF_NFLG'):
        # skip consensus sequences
        continue

    patid, coldate, sample, snum, stype = qname.split('_')
    dt = datetime.strptime(coldate, "%d%b%Y")

    if patid not in by_patient:
        by_patient.update({patid: []})
    by_patient[patid].append({'accno': accno, 'coldate': dt, 'sample': sample, 'seq': s})


for patid, data in by_patient.items():
    outfile = open('data/{}.fa'.format(patid), 'w')
    for row in data:
        if row['sample'] in ['SC', '1Y', 'CP']:
            moltype = 'RNA'
        elif row['sample'] in ['CC', '1A', '2A']:
            moltype = 'DNA'
        else:
            print("ERROR: unrecognized sample type {}".format(row['sample']))
            sys.exit()
        h = "{}_{}_{}_{}".format(row['accno'], row['sample'], moltype,
                              row['coldate'].date().isoformat())
        outfile.write(">{}\n{}\n".format(h, row['seq']))
    outfile.close()
