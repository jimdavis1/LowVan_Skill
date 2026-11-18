#!/usr/bin/env python3
"""Give each bar the (label, value, pct, note) shape gen_coverage_audit expects,
and attach the count behind each percentage."""
import json, os
p = os.path.expanduser('~/lowvan_scratch/tobani/pages/tobaniviridae_audit.json')
d = json.load(open(p))
NOTES = {
 "Torovirinae S / M / N / ORF1a / Mpro": "24 of 24 panel genomes",
 "Torovirinae HE": "23 of 24; the miss is a genome whose HE region is truncated",
 "Torovirinae ORF1ab": "21 of 24; every miss carries ambiguous bases",
 "Piscanivirinae S / M / N / Mpro": "13 of 13",
 "Piscanivirinae ORF1a": "11 of 13; two contigs are 5'-truncated past pp1a",
 "Piscanivirinae ORF1ab": "8 of 13; see section 07",
 "Remotovirinae S / M / N / ORF1a / ORF6": "8 of 8",
 "Remotovirinae ORF1ab": "7 of 8; the miss carries 14 ambiguous bases",
 "Serpentovirinae S / M": "14 of 14",
 "Serpentovirinae N / ORF1a": "13 of 14",
 "Serpentovirinae ORF1ab": "12 of 14",
}
for s in d['sections']:
    for b in s['blocks']:
        if b['type'] == 'bars':
            b['unit'] = '%'
            b['rows'] = [[lab, '%g' % pct, pct, NOTES.get(lab, '')]
                         for lab, pct in b['rows']]
json.dump(d, open(p, 'w'), indent=1)
print('bars fixed')
