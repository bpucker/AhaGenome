### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """
					python3 nummeric_check.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					bug reports and feature requests: b.pucker@tu-braunschweig.de
					"""

import os, sys, re, glob, subprocess
from operator import itemgetter

# --- end of imports --- #


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	min_marker_cutoff = 2
	gen_diff_cutoff = 0.1
	
	# --- load marker from file --- #
	markers = []
	marker_order = {}
	counter = 0
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			counter += 1
			parts = line.strip().split('\t')
			markers.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'line': line, 'genpos': float( parts[3].split(':')[-1] ) } )
			marker_order.update( { line: counter } )
			line = f.readline()
	
	potential_errors = []
	for idx, marker in enumerate( markers ):
		if idx > 0:
			if marker['chr'] == markers[ idx-1 ]['chr']:
				if markers[ idx-1 ]['pos'] > marker['pos']:
					potential_errors.append( marker )
	
	print ( len( potential_errors ) )
	
	# --- write marker into output file --- #
	block_counter = 0
	with open( output_file, "w" ) as out:
		index = 0
		while True:
			tmp_errors = []
			while marker_order[ potential_errors[ index ]['line'] ]+1 == marker_order[ potential_errors[ index+1 ]['line'] ]:
				tmp_errors.append( potential_errors[ index ] )
				index += 1
				if index >= len( potential_errors )-1:
					break
			if len( tmp_errors ) > min_marker_cutoff:
				gen_diff = abs( tmp_errors[0]['genpos'] - tmp_errors[-1]['genpos'] )
				if gen_diff > gen_diff_cutoff:
					for marker in tmp_errors:
						out.write( marker['line'] )
					out.write( "\n" )
					block_counter += 1
			index += 1
			if index >= len( potential_errors )-1:
				break
	
	print ( block_counter )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
