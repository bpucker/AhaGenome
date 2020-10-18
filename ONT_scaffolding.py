### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###


__usage__ = """
					python ONT_scaffolding.py
					--assembly <ASSEMBLY_FILE>
					--reads <READS_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--minlen <MINIMAL_ALIGNMENT_LENGTH>[500]
					--minsim <MINIMAL_ALIGNMENT_SIMILARITY>[80]
					--cpus <THREADS_FOR_BLAST_SEARCH>[20]
					"""

import os, sys

# --- end of imports --- #


def load_BLAST_results( result_file, len_cutoff, sim_cutoff ):
	"""! @brief load best BLAST hit per subject """
	
	#reads are keys in outer dictionary
	#contig names are keys in inner dictionaries
	BLAST_hits = {}
	with open( result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if int( parts[3] ) > len_cutoff:
				if float( parts[2] ) > sim_cutoff:
					try:
						if float( parts[-1] ) > BLAST_hits[ parts[1] ][ parts[0] ]['score']:
							BLAST_hits[ parts[1] ][ parts[0] ] = { 'pos': ( int( parts[8] ) + int( parts[9] ) ) / 2.0, 'score': float( parts[-1] ) }
					except KeyError:
						try:
							BLAST_hits[ parts[1] ].update( { parts[0]: { 'pos': ( int( parts[8] ) + int( parts[9] ) ) / 2.0, 'score': float( parts[-1] ) } } )
						except KeyError:
							BLAST_hits.update( { parts[1]: { parts[0]: { 'pos': ( int( parts[8] ) + int( parts[9] ) ) / 2.0, 'score': float( parts[-1] ) } } } )
			line = f.readline()
	return BLAST_hits


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


def main( parameters ):
	"""! @brief run everything """
	
	assembly_file = parameters[ parameters.index("--assembly")+1 ]
	read_file = parameters[ parameters.index("--reads")+1 ]
	output_folder = parameters[ parameters.index("--out")+1 ]
	
	if "--minlen" in parameters:
		len_cutoff = int( parameters[ parameters.index("--minlen")+1 ] )
	else:
		len_cutoff = 500
	
	if '--minsim' in parameters:
		sim_cutoff = int( parameters[ parameters.index("--minsim")+1 ] )
	else:
		sim_cutoff = 80
	
	if '--cpus' in parameters:
		cpus = int( parameters[ parameters.index("--cpus")+1 ] )
	else:
		cpus = 20
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	
	# --- prepare conig end file --- #
	contig_end_file = output_folder + "contig_ends.fasta"
	if not os.path.isfile( contig_end_file ):
		seqs = load_sequences( assembly_file )
		with open( contig_end_file, "w" ) as out:
			for key in seqs.keys():
				out.write( '>' + key + "_1\n" + seqs[ key ][:10000] + "\n" )
				out.write( '>' + key + "_2\n" + seqs[ key ][-10000:] + "\n" )
	
	# --- construct database from reads --- #
	db = output_folder + "blast_db"
	if not os.path.isfile( db + ".nal" ):
		os.popen( "makeblastdb -in " + read_file + " -out " + db + " -dbtype nucl" )
	
	# --- BLASTn contig ends vs. db --- #
	blast_result_file = output_folder + "BLAST_results.txt"
	if not os.path.isfile( blast_result_file ):
		os.popen( "blastn -query " + contig_end_file + " -db " + db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.001 -num_threads " + str( cpus ) )
	
	# --- process BLAST results --- #
	blast_hits = load_BLAST_results( blast_result_file, len_cutoff, sim_cutoff )
	
	info_file = output_folder + "info.txt"
	matches = []
	with open( info_file, "w" ) as out:
		for read in blast_hits.keys():
			data = blast_hits[ read ]
			contigs = data.keys()
			matches.append( contigs )
			out.write( read + "\t" + ",".join( contigs ) + "\n" )
	
	summary_file = output_folder + "summary.txt"
	with open( summary_file, "w" ) as out:
		contig_end_IDs = load_sequences( contig_end_file ).keys()
		for idx1, c1 in enumerate( contig_end_IDs ):
			for idx2, c2 in enumerate( contig_end_IDs ):
				if idx2 > idx1:
					counter = 0
					for match in matches:
						if c1 in match and c2 in match:
							counter += 1
					out.write( c1 + "\t" + c2 + "\t" + str( counter ) + "\n" )


if "--assembly" in sys.argv and '--reads' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
