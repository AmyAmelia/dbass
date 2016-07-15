"""
Label events as de novo (i.e. exon extension/truncation), cryptic, pseudoexon, or unclear.
"""

__author__ = 'bernie'

import re
import sys
import argparse


def get_authentic_donor_idx(seq, new_donor_idx, 
							donor_regex="[AGCT][agct]",
							acceptor_regex="[acgt][ACGT]",
							get_match_object=False):
	authentic_donors = re.finditer(donor_regex, seq)
	authentic_acceptors = re.finditer(acceptor_regex, seq)
	prev_ag = None
	next_ag = None
	the_donor = None
	for acceptor in authentic_acceptors:
		if new_donor_idx < acceptor.start():
			next_ag = acceptor.start()
			break
		prev_ag = acceptor.start()
	for donor in authentic_donors:
		if donor.start() > prev_ag:
			the_donor = donor if get_match_object else donor.start()
			break
	return the_donor, prev_ag, next_ag


def get_authentic_acceptor_idx(seq, new_acceptor_idx, 
							   donor_regex="[AGCT][agct]",
							   acceptor_regex="[acgt][ACGT]",
							   get_match_object=False):
	authentic_donors = re.finditer(donor_regex, seq)
	authentic_acceptors = re.finditer(acceptor_regex, seq)
	prev_gt = None
	next_gt = None
	the_acceptor = None
	for donor in authentic_donors:
		if new_acceptor_idx < donor.start():
			next_gt = donor.start()
			break
		prev_gt = donor.start()
	for acceptor in authentic_acceptors:
		if acceptor.start() > prev_gt:
			the_acceptor = acceptor if get_match_object else acceptor.start()
			break
	return the_acceptor, prev_gt, next_gt


def label_donor(row, header):
	if re.search('[Pp]seudo ?exon', row[header['Comment']]):
		return 'pseudoexon'
	
	seq = row[header['NucleotideSequence']]
	event_regex = "[\[\(].*[\]\)]"
	if len(re.findall('/', seq)) != 1 or len(re.findall(event_regex, seq)) != 1:
		return 'unclear'
	else:
		new = re.search('/', seq).start()
		event = re.search(event_regex, seq)

	# find the authentic donor site that gets outcompted by an aberrant site
	donor_regex = "[A-Z]([\[\(].*[\]\)])?[a-z]"
	acceptor_regex = "[a-z]([\[\(].*[\]\)])?[A-Z]"
	results = get_authentic_donor_idx(seq, new, donor_regex=donor_regex, acceptor_regex=acceptor_regex, get_match_object=True)
	donor, exon_start, intron_end = results
	
	if not donor: return 'insufficient'
	
	alteration = row[header['Alteration']]
	donor_str = seq[donor.start():donor.end()]
	event_str = seq[event.start():event.end()]
	if event_str in donor_str: return 'cryptic'
	
	if event.start() > donor.end():
		between = seq[donor.end():event.start()]
		dist = len(filter(lambda c:re.match('[ATGCatgc]', c), between)) + 2
	else:
		between = seq[event.end():donor.start()]
		dist = -len(filter(lambda c:re.match('[ATGCatgc]', c), between)) - 2 
	if dist >= -3 and dist <= 6: return 'cryptic'
	
	if event.start() >= new:
		between = seq[new:event.start()]
		dist = len(filter(lambda c:re.match('[ATGCatgc]', c), between))
	else:
		between = seq[event.end():new]
		dist = -len(filter(lambda c:re.match('[ATGCatgc]', c), between))
	if dist < -2 or dist > 5: return 'trans'
	
	return 'denovo'


def main(args):
	o = open(args.output, 'w') if args.output != sys.stdout else sys.stdout
	o_header = ['Gene', 'Alteration', 'NucleotideSequence', 'Label', 'SpliceSiteType']
	o.write('\t'.join(o_header) + '\n')
	header = None
	i = 0
	dbass = sys.stdin if args.input == '-' else open(args.input, 'r') 
	for row in dbass:
		if not header:
			fields = row.strip().split('\t')
			header = dict(zip(fields, range(len(fields))))
			continue
		row = row.strip().split('\t')
		label = label_donor(row, header)
		newrow = [row[header['GeneName']], row[header['Alteration']], row[header['NucleotideSequence']], label, args.splice]
		o.write('\t'.join(newrow) + '\n')
	o.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', '--in', '-i', dest='input', default='-', help='Input DBASS file. Accepts stdin as default')
	parser.add_argument('--output', '--out', '-o', dest='output', default=sys.stdout)
	parser.add_argument('-s', dest='splice', required=True, choices=['acceptor', 'donor'], default='donor')
	args = parser.parse_args()

	main(args)

