#!/usr/bin/python2

from Bio import SeqIO

gb=SeqIO.parse('pedv-usa-mexico.gb','genbank')

recs = []
for g in gb:
    for feat in g.features:
        if feat.type == 'source':
            if 'collection_date' in feat.qualifiers and 'country' in feat.qualifiers:
                rec =  [str(g.id), feat.qualifiers['country'][0]]
                rec.append(feat.qualifiers['collection_date'][0])
                rec = ','.join(rec)
                recs.append(rec)
            
f = open('pedv-usa-mexico.csv', 'w')
f.write('accession,country,date\n')
for rec in recs:
    f.write(rec + '\n')
f.close()
