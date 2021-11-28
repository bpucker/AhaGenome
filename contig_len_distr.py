### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python contig_len_distr.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					"""

from operator import itemgetter
import sys, os

# --- end of imports --- #

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	assembly = load_sequences( input_file )
	data = []
	total = 0.0
	for key in assembly.keys():
		data.append( { 'ID': key, 'len': len( assembly[ key ] ) } )
		total += len( assembly[ key ] )

	counter = 0
	with open( output_file, "w" ) as out:
		out.write( "#\tContig\tContigLength\tCumulativeLength\tContribution[%]\tCumulativeContribution[%]\n" )
		for idx, each in enumerate( sorted( data, key=itemgetter('len') )[::-1] ):
			counter += each['len']
			out.write( "\t".join( map( str, [ idx+1, each['ID'], each['len'], counter, 100*each['len']/total, 100*counter/total ] ) ) + "\n" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
