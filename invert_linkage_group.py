### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """
					python3 invert_linkage_group.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					"""

from operator import itemgetter
import os, sys

# --- end of imports --- #

def load_genetic_map( input_file ):
	"""! @brief load genetic map from BED file """
	
	genetic_map = {}
	lg_order = []
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				genetic_map[ parts[3].split(':')[0] ].append( { 'contig': parts[0], 'start': parts[1], 'end': parts[2], 'genpos': float( parts[3].split(':')[1] ), 'contigpos': parts[4] } )
			except KeyError:
				genetic_map.update( { parts[3].split(':')[0]: [ { 'contig': parts[0], 'start': parts[1], 'end': parts[2], 'genpos': float( parts[3].split(':')[1] ), 'contigpos': parts[4] } ] } )
				lg_order.append( parts[3].split(':')[0] )
			line = f.readline()
	sorted_genetic_map = {}
	for lg in lg_order:
		sorted_genetic_map.update( { lg: sorted( genetic_map[ lg ], key=itemgetter('genpos') ) } )
	return sorted_genetic_map, lg_order


def calculate_inverted_linkage_groups( genetic_map ):
	"""! @brief calculate an inverted orientation of genetic markers per linkage group """
	
	inverted_genetic_map = {}
	for lg in list( genetic_map.keys() ):
		# --- load distances of markers in linkage map and on contig --- #
		markers_with_distance = []
		for idx, marker in enumerate( genetic_map[ lg ] ):
			try:
				marker.update( { 'dist': genetic_map[ lg ][ idx+1 ]['genpos'] - marker['genpos'] } )
			except IndexError:
				marker.update( { 'dist': 0 } )
			markers_with_distance.append( marker )
		
		# ---- calculate new genpos and contigpos --- #
		current_lg_pos = 0.0
		current_contig_pos = 0.0
		fin_marker_data = []
		for marker in markers_with_distance[::-1]:	#go through inverted list and calculate positions
			marker['genpos'] = current_lg_pos + marker['dist']
			fin_marker_data.append( marker )
			current_lg_pos += marker['dist']
		inverted_genetic_map.update( { lg: fin_marker_data } )
	return inverted_genetic_map


def main( arguments ):
	"""! @brief handle everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]

	genetic_map, lg_order = load_genetic_map( input_file )
	inverted_linkage_group_map = calculate_inverted_linkage_groups( genetic_map )


	with open( output_file, "w" ) as out:
		for lg in lg_order:
			markers = sorted( inverted_linkage_group_map[ lg ], key=itemgetter('genpos') )
			for marker in markers:
				out.write( "\t".join( [	marker['contig'],
													marker['start'],
													marker['end'],
													lg + ":" + str( marker['genpos'] ),
													marker['contigpos']
													] ) + "\n" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
