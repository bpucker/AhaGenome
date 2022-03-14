### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """
					python3 manual_scaffolding.py
					--fasta <FASTA_INPUT_FILE>
					--agp <AGP_INPUT_FILE>
					--out <OUTPUT_FASTA_FILE>
					"""

from operator import itemgetter
import os, sys

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


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def load_agp_file( agp_input_file, assembly ):
	"""! @brief load information from AGP file """
	
	chr_structures = {}
	chr_order = []
	with open( agp_input_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if parts[4] == "W":
					seq = assembly[ parts[5] ]
					if parts[-1] == "-":
						seq = revcomp( seq )
					try:
						chr_structures[ parts[0] ].append( seq )
					except KeyError:
						chr_structures.update( { parts[0]: [ seq ] } )
						chr_order.append( parts[0] )
				else:	#add gap
					chr_structures[ parts[0] ].append( "N"*int( parts[5] ) )	#ADD THIS
			line = f.readline()
	return chr_structures, chr_order


def main( arguments ):
	"""! @brief handle everything """
	
	fasta_input_file = arguments[ arguments.index('--fasta')+1 ]
	agp_input_file = arguments[ arguments.index('--agp')+1 ]
	output_fasta_file = arguments[ arguments.index('--out')+1 ]

	assembly = load_sequences( fasta_input_file )
	chromosome_structure, lg_order = load_agp_file( agp_input_file, assembly )

	with open( output_fasta_file, "w" ) as out:
		for lg in lg_order:
			seqs = chromosome_structure[ lg ]
			out.write( ">" + lg + "\n" + "".join( seqs ) + "\n" )


if '--fasta' in sys.argv and '--agp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
