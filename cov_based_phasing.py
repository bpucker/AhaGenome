### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """
					python3 cov_based_phasing.py
					--covinfo1 <COVERAGE_INFO_FILE1>
					--covinfo2 <COVERAGE_INFO_FILE1>
					--assembly <ASSEMBLY_FILE>
					--haplo1 <HAPLO1_OUTPUT_FILE>
					--haplo2 <HAPLO2_OUTPUT_FILE>
					--unclass <UNCLASSIFIED_OUTPUT_FILE>
					bug reports and feature requests: b.pucker@tu-braunschweig.de
					"""

import os, sys

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


def load_avg_cov( cov_info_file_D111 ):
	"""! @brief load average coverage value per contig """
	
	infos = {}
	with open( cov_info_file_D111, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			infos.update( { parts[0]: float( parts[2] ) } )
			line = f.readline()
	
	return infos
	
	
def main( arguments ):
	"""! @brief run everything """

	cov_info_file_D111 = arguments[ arguments.index('--covinfo1')+1 ]
	cov_info_file_D654 = arguments[ arguments.index('--covinfo2')+1 ]
	assembly_file = arguments[ arguments.index('--assembly')+1 ]

	D111_output_file = arguments[ arguments.index('--haplo1')+1 ]
	D654_output_file = arguments[ arguments.index('--haplo2')+1 ]
	unclass_output_file = arguments[ arguments.index('--unclass')+1 ]

	factor = 10

	assembly = load_sequences( assembly_file )

	cov_values_D111 = load_avg_cov( cov_info_file_D111 )
	cov_values_D654 = load_avg_cov( cov_info_file_D654 )


	D111 = 0
	D654 = 0
	unclass = 0

	with open( D111_output_file, "w" ) as out1:
		with open( D654_output_file, "w" ) as out2:
			with open( unclass_output_file, "w" ) as out3:
				for key in assembly.keys():
					D111_cov = cov_values_D111[ key ]
					D654_cov = cov_values_D654[ key ]
					if D111_cov > factor*D654_cov:
						D111 += len( assembly[ key ] )
						out1.write( '>' + key + "\n" + assembly[ key ] + "\n" )
					elif D654_cov > factor*D111_cov:
						D654 += len( assembly[ key ] )
						out2.write( '>' + key + "\n" + assembly[ key ] + "\n" )
					else:
						print ( key + "\t" + str( len( assembly[ key ] ) ) + "\t" + str( D111_cov ) + "\t" + str( D654_cov ) )
						unclass += len( assembly[ key ] )
						out3.write( '>' + key + "\n" + assembly[ key ] + "\n" )
					

			print ( "haplo1: " + str( D111 ) )	#D111
			print ( "haplo2: " + str( D654 ) )	#D654
			print ( "unclass: " + str( unclass ) )


if '--covinfo1' in sys.argv and '--covinfo2' in sys.argv and '--assembly' in sys.argv and '--haplo1' in sys.argv and '--haplo2' in sys.argv and '--unclass' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
	
