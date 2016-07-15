"""
Mutate NucleotideSequence field of DBASS data to reflect alternate allele.
"""

import sys
import re
import argparse
import fileinput


def main(args):
	o = open(args.output, 'w') if args.output != sys.stdout else sys.stdout
	i = 1
	dbass = sys.stdin if args.input == '-' else open(args.input, 'r')
	for row in dbass:
		fields = row.strip().split('\t')
		if i == 1:
			header = dict(zip(fields, range(len(fields))))
			o.write('\t'.join(fields) + '\n')
			i += 1
		else:
			orig = fields[header['NucleotideSequence']]
			seq = re.sub('[\[\]]', '', orig) # include insertions
			delta = 0
			for event in re.finditer('\(.*?\)', seq):
				l = len(seq)
				event_str = event.group(0)
				start = event.start() + delta
				end = event.end() + delta
				# SNP
				if '>' in event_str:
					match = re.search('\((.*)>(.*)\)', event_str)
					allele = match.group(1) if args.ref else match.group(2) 
					seq = seq[:start] + allele + seq[end:]
				# deletion
				elif args.ref:
					seq = filter(lambda char: char not in ['(', ')'], seq)
				else:
					seq = seq[:start] + seq[end:]
				delta = len(seq) - l

			newfields = fields
			newfields[header['NucleotideSequence']] = seq
			o.write('\t'.join(newfields) + '\n')	

	o.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--dbass', '-i', dest='input', default='-', help='Input DBASS file. Accepts stdin as default.')
	parser.add_argument('--out', '-o', dest='output', default=sys.stdout)
	parser.add_argument('--ref', action='store_true', help='Derive wild-type sequence instead')
	args = parser.parse_args()
	main(args)

		

