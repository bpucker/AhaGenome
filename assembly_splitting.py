### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.2 ###

__usage__ = """
					python3 assembly_splitting.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					--info <INFO_FILE>
					bug reports and feature requests: b.pucker@tu-braunschweig.de
					"""

import os, sys
from operator import itemgetter

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(" ")[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(" ")[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_infos( info_file ):
	"""! @brief load infos from file """
	
	infos = {}
	with open( info_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				infos[ parts[1] ].append( { 'start': int( parts[2] ), 'end': int( parts[3] ) } )
			except KeyError:
				infos.update( { parts[1]: [ { 'start': int( parts[2] ), 'end': int( parts[3] ) } ] } )
			line = f.readline()
	return infos


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	info_file = arguments[ arguments.index('--info')+1 ]
	
	min_seq_len = 10000
	
	assembly = load_sequences( input_file )
	info = load_infos( info_file )
	
	with open( output_file, "w" ) as out:
		for key in list( assembly.keys() ):
			counter = 1
			try:
				todo = sorted( info[ key ], key=itemgetter( 'start' ) )
			except KeyError:
				todo = []
				out.write( ">" + key + "\n" + assembly[ key ] + "\n" )
			if len( todo ) > 0:
				if len( todo ) == 1:
					seq = assembly[ key ][ :todo[0]['start'] ]
					if len( seq ) > min_seq_len:
						out.write( ">" + key + str( counter ) + "\n" + seq + "\n" )
						counter += 1
					seq = assembly[ key ][ todo[0]['end']: ]
					if len( seq ) > min_seq_len:
						out.write( ">" + key + str( counter ) + "\n" + seq + "\n" )
						counter += 1
				else:
					for idx, each in enumerate( todo ):
						if idx == 0:
							seq = assembly[ key ][ :todo[0]['start'] ]
							if len( seq ) > min_seq_len:
								out.write( ">" + key + str( counter ) + "\n" + seq + "\n" )
								counter += 1
						if idx > 0:
							seq = assembly[ key ][ todo[idx-1]['end']:todo[idx]['start'] ]
							if len( seq ) > min_seq_len:
								out.write( ">" + key + str( counter ) + "\n" + seq + "\n" )
								counter += 1
						if idx == len( todo )-1:
							seq = assembly[ key ][ todo[idx]['end']: ]
							if len( seq ) > min_seq_len:
								out.write( ">" + key + str( counter ) + "\n" + seq + "\n" )
								counter += 1
			


if '--in' in sys.argv and '--out' in sys.argv and '--info' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
