### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python genetic_map_constructor.py
					--csv <GENETIC_MAP_IN_CSV>
					--ref <REFERENCE_ASSEMBLY>
					--assembly <ASSEMBLY_FOR_SCAFFOLDING>
					--out <OUTPUT_FOLDER>
					"""


import sys, os
from operator import itemgetter

# --- end of imports --- #



def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_genetic_map( genetic_map_file ):
	"""! @brief load genetic map from given file """
	
	with open( genetic_map_file, "r" ) as f:
		genomic_pos = f.readline().strip().split(',')[1:]
		linkage_group = f.readline().strip().split(',')[1:]
		genetic_pos = f.readline().strip().split(',')[1:]
	
	genetic_map = {}
	for idx, pos in enumerate( genomic_pos ):
		genetic_map.update( { pos: { 'group': linkage_group[ idx ], 'pos': genetic_pos[ idx ] } } )	#genetic map position per genomic position
	return genetic_map


def load_best_hit( blast_result_file ):
	"""! @brief load best BLAST hit per query """
	
	BLAST_hits  ={}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if float( parts[-1] ) > BLAST_hits[ parts[0] ]['score']:
					BLAST_hits[ parts[ 0 ] ] = { 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )/ 2.0, 'score': float( parts[-1] ) }
			except KeyError:
				BLAST_hits.update( { parts[0]: { 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )/ 2.0, 'score': float( parts[-1] ) } } )
			line = f.readline()
	return BLAST_hits


def construct_bed_file( bed_file_for_allpaths, best_hit_per_query, genetic_map ):
	"""! @brief construct BED file with genetic markers """
	
	data = []
	for each in genetic_map.keys():
		try:
			assembly_pos = best_hit_per_query[ each ]
			data.append( { 	'chr': best_hit_per_query[each]['chr'],
										'pos1': str( int( best_hit_per_query[each]['pos']-1 ) ),
										'pos2': str( int( best_hit_per_query[each]['pos'] ) ),
										'genetic_pos': "X-" + str( genetic_map[each]['group'] ) + ":" + str( genetic_map[each]['pos'] ),
										'genomic_pos': str( best_hit_per_query[each]['chr'] ) + ":" + str( best_hit_per_query[each]['pos'] ),
										'group': genetic_map[each]['group'],
										'pos': genetic_map[each]['pos']
									} )
		except KeyError:
			pass
	sorted_data = sorted( data, key=itemgetter('group', 'pos') )
	with open( bed_file_for_allpaths, "w" ) as out:
		for point in sorted_data:
			new_line = [ 	point['chr'],
									point['pos1'],
									point['pos2'],
									point['genetic_pos'],
									point['genomic_pos']
								]
			out.write( "\t".join( map( str, new_line ) ) + "\n" )


def main( arguments ):
	"""! @brief handling files and data """


	genetic_map_file = arguments[ arguments.index('--csv')+1 ]
	genetic_map_reference = arguments[ arguments.index('--ref')+1 ]
	
	assembly_for_scaffolding = arguments[ arguments.index('--assembly')+1 ]
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	genetic_map = load_genetic_map( genetic_map_file )
	ref_assembly = load_sequences( genetic_map_reference )
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- extract 500 bp left +right --- #
	query_file = output_folder + "query.fasta"
	if not os.path.isfile( query_file ):
		with open( query_file, "w" ) as out:
			for marker in genetic_map.keys():
				seq = ref_assembly[ marker.split('_')[0] ][ int( marker.split('_')[1] )-500:int( marker.split('_')[1] )+500 ]
				out.write( '>' + marker + "\n" + seq + "\n" )
		
	# --- run BLAST --- #
	blast_result_file = output_folder + "blast_results.txt"
	if not os.path.isfile( blast_result_file ):
		os.popen( "blastn -query " + query_file + " -subject " + assembly_for_scaffolding + " -out " + blast_result_file + " -outfmt 6 -evalue 0.0001"  )
	
	
	# --- find best hit --- #
	best_hit_per_query = load_best_hit( blast_result_file )
	
	# --- construct BED file --- #
	bed_file_for_allpaths = output_folder + "genetic_markers.bed"
	construct_bed_file( bed_file_for_allpaths, best_hit_per_query, genetic_map )


if '--csv' in sys.argv and '--ref' in sys.argv and '--assembly' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

