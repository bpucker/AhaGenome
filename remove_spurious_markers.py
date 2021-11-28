### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """
					python3 remove_spurious_markers.py
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
	
	cutoff = 1000000
	gen_dist_cutoff = 0.0001
	
	# --- load marker from file --- #
	markers = []
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			markers.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'line': line, 'genpos': float( parts[3].split(':')[-1] ) } )
			line = f.readline()
	
	clean_markers = []
	for idx, marker in enumerate( markers ):
		if idx > 0 and idx < len( markers )-1:
			status = True
			if marker['chr'] != markers[ idx-1 ]['chr']:
				if marker['chr'] != markers[ idx+1 ]['chr']:
					status = False
				elif marker['pos'] > markers[ idx-1 ]['pos']+cutoff:
					if marker['pos'] < markers[ idx+1 ]['pos']-cutoff:
						status = False
			if marker['pos'] > markers[ idx-1 ]['pos']+cutoff:
				if marker['pos'] < markers[ idx+1 ]['pos']-cutoff:
					status = False
			if status:
				if len( clean_markers ) > 0:
					if clean_markers[-1]['chr'] == marker['chr']:
						if abs( clean_markers[-1]['genpos'] - marker['genpos'] ) >= gen_dist_cutoff:
							clean_markers.append( marker )
					else:
						clean_markers.append( marker )
				else:
					clean_markers.append( marker )
	
	
	print ( len( markers ) )
	print ( len( clean_markers ) )
	
	# --- write marker into output file --- #
	with open( output_file, "w" ) as out:
		for marker in clean_markers:
			out.write( marker['line'] )




if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
