### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python phasing_check.py
					--ref <ASSEMBLY_FILE>
					--vcf1 <VCF_FILE1>
					--vcf2 <VCF_FILE2>
					--out <OUTPUT_FOLDER>
					
					optional:
					--covinfo1 <COV_PER_CONTIG_FILE1>
					--covinfo2 <COV_PER_CONTIG_FILE2>
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


def load_variants( vcf, min_cov=10, min_freq=0.1 ):
	"""! @brief load variants """
	
	variants = {}
	frequencies = []
	coverages = []
	with open( vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if len( parts[3] ) == 1 and len( parts[4] ) == 1:
					if parts[-1] != "./.":
						allele_freqs = map( int, parts[-1].split(':')[1].split(',') )
						if len( allele_freqs ) == 2:
							if sum( allele_freqs ) > min_cov:	#valid variant
								frequencies.append( min( allele_freqs ) / float( sum( allele_freqs ) ) )
								coverages.append( sum( allele_freqs ) )
								if min( allele_freqs ) / float( sum( allele_freqs ) ) > min_freq:
									status = "het"
								else:
									status = "homo"
								try:
									variants[ parts[0] ].append( { 'pos': int( parts[1] ), 'status': status } )
								except KeyError:
									variants.update( { parts[0]: [ { 'pos': int( parts[1] ), 'status': status } ] } )
			line = f.readline()
	return variants, frequencies, coverages


def count_homo_and_hetero_vars( variants ):
	"""! @brief count homozygous and heterozygous variants in given list """
	
	homo = 0
	hetero = 0
	for each in variants:
		if each['status'] == "het":
			hetero += 1
		else:
			homo += 1
	return homo, hetero


def analyze_data( frequencies1, frequencies2, coverages1, coverages2, output_folder ):
	"""! @brief analyse data sets """
	
	import matplotlib.pyplot as plt
	
	freq_fig_file = output_folder + "freq.png"
	fig, ax = plt.subplots()
	ax.hist( frequencies1, bins=100, color="lime", label="freq1", alpha=0.5 )
	ax.hist( frequencies2, bins=100, color="magenta", label="freq2", alpha=0.5 )
	ax.set_xlim( 0, 1 )
	ax.set_xlabel( "minor variant frequency" )
	ax.set_ylabel( "number of variants" )
	fig.savefig( freq_fig_file, dpi=300 )
	
	cov_fig_file = output_folder + "cov.png"
	fig, ax = plt.subplots()
	ax.hist( coverages1, bins=10000, color="lime", label="cov1", alpha=0.5 )
	ax.hist( coverages2, bins=10000, color="magenta", label="cov2", alpha=0.5 )
	ax.set_xlim( 0, 50 )
	ax.set_xlabel( "coverage" )
	ax.set_ylabel( "number of variants" )
	fig.savefig( cov_fig_file, dpi=300 )


def load_cov_per_contig( cov_info_file ):
	"""! @brief load average coverage value per contig """
	
	cov_per_contig = {}
	with open( cov_info_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			cov_per_contig.update( { parts[0]: parts[2] } )
			line = f.readline()
	return cov_per_contig


def generate_var_distr_plots( variants1, variants2, assembly, min_contig_len, window_size, output_folder ):
	"""! @brief generate a plot per contig """
	
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches
	
	for key in assembly.keys():
		if len( assembly[ key ] ) > min_contig_len:
			fig_file = output_folder + key + ".png"
			
			# --- prepare data structure --- #
			x = []
			y_homo1 = []
			y_hetero1 = []
			y_hetero2 = []
			y_homo2 = []
			
			for i in range( ( len( assembly[ key ] ) / window_size ) +1 ):
				x.append( ( ( i+0.5 ) * window_size ) / 1000000.0 )
				y_homo1.append( 0 )
				y_homo2.append( 0 )
				y_hetero1.append( 0 )
				y_hetero2.append( 0 )
			
			# --- count variants per window --- #
			try:
				for var in variants1[ key ]:
					if var['status'] == "homo":
						y_homo1[ var['pos'] / window_size ] += 1
					else:
						y_hetero1[ var['pos'] / window_size ] += 1
			except KeyError:
				pass
			
			try:
				for var in variants2[ key ]:
					if var['status'] == "homo":
						y_homo2[ var['pos'] / window_size ] += 1
					else:
						y_hetero2[ var['pos'] / window_size ] += 1
			except KeyError:
				pass
			
			ratios1 = []
			ratios2 = []
			counts1 = []
			for idx, val in enumerate( y_homo1 ):
				ratios1.append( float( val ) / ( 1 + val +y_hetero1[ idx ]  ) )
				counts1.append( y_homo1[ idx ] + y_hetero1[ idx ] )
				ratios2.append( float( y_homo2[ idx ] ) / ( 1 + y_homo2[ idx ] +y_hetero2[ idx ]  ) )
			
			# --- generate figure -- #
			fig, ax = plt.subplots( figsize=( 20, 4 ) )
			
			ax.plot( 0, len( assembly[ key ] ), color="black" )
			
			ax.scatter( x, ratios1, color="magenta", s=1 )
			ax2 = ax.twinx()
			ax2.scatter( x, ratios2, color="lime", s=1 )
			ax3 = ax2.twinx()
			ax3.scatter( x, counts1, s=1, color="black" )
			
			ax.set_xlim( 0, max( x ) )
			
			ax.set_ylim( 0, 1 )
			ax2.set_ylim( 0, 1 )
			ax3.set_ylim( 0, max( counts1 ) )
			
			legend = [	mpatches.Patch(color='magenta', label='D111'),	#add command line argument to specify names
								mpatches.Patch(color='lime', label='D654'),
								mpatches.Patch(color='black', label='D111 counts')
							]
			ax.legend( handles=legend, fontsize=10 )
			
			# ax.scatter( x, y_homo1, s=3, label="homo1", color="orange", alpha=0.5 )
			# ax2 = ax.twinx()
			# ax2.scatter( x, y_hetero1, s=3, label="hetero1", color="red", alpha=0.5 )
			
			# ax3 = ax.twinx()
			# ax3.scatter( x, y_homo2, s=3, label="homo2", color="lime", alpha=0.5 )
			# ax4 = ax.twinx()
			# ax4.scatter( x, y_hetero2, s=3, label="hetero2", color="green", alpha=0.5 )
			
			# legend = [	mpatches.Patch(color='orange', label='homo1'),
								# mpatches.Patch(color='red', label='hetero1'),
								# mpatches.Patch(color='lime', label='homo2'),
								# mpatches.Patch(color='green', label='hetero2')
							# ]
			
			# ax.legend( handles=legend, fontsize=10 )
			
			# ax.set_ylim( 0, max( y_homo1 ) )
			
			ax.set_ylabel( "proportion of homozygot variants" )
			ax.set_xlabel( "position on contig [Mbp]" )
			
			plt.subplots_adjust( left=0.03, right=0.98, top=0.999, bottom=0.15 )
			
			fig.savefig( fig_file, dpi=300 )
			plt.close( "all" )


def main( arguments ):
	"""! @brief run everything """
	
	assembly_file = arguments[ arguments.index('--ref')+1 ]
	vcf1 = arguments[ arguments.index('--vcf1')+1 ]
	vcf2 = arguments[ arguments.index('--vcf2')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	output_file = output_folder + "SUMMARY.txt"
	
	if '--size' in arguments:
		try:
			min_contig_len = int( arguments[ arguments.index('--size')+1 ] )
		except:
			min_contig_len = 500000
	else:
		min_contig_len = 500000

	if '--covinfo1' in arguments:
		cov_info_file1 = arguments[ arguments.index('--covinfo1')+1 ]
		cov_info1 = load_cov_per_contig( cov_info_file1 )
	else:
		cov_info1 = {}

	if '--covinfo2' in arguments:
		cov_info_file2 = arguments[ arguments.index('--covinfo2')+1 ]
		cov_info2 = load_cov_per_contig( cov_info_file2 )
	else:
		cov_info2 = {}
	
	if '--window' in arguments:
		window_size = int( arguments[ arguments.index('--window')+1 ] )
	else:
		window_size = 100000
	
	
	# --- loading all data sets --- #
	variants1, frequencies1, coverages1 = load_variants( vcf1, min_cov=10, min_freq=0.1 )
	variants2, frequencies2, coverages2 = load_variants( vcf2, min_cov=10, min_freq=0.1 )
	
	assembly = load_sequences( assembly_file )
	
	# --- generating summary file --- #
	with open( output_file, "w" ) as out:
		out.write( "contig\thomo1\thetero1\tcoverage1\thomo2\thetero2\tcoverage2\tcontig_length\n" )
		for key in assembly.keys():
			if len( assembly[ key ] ) > min_contig_len:
				try:
					cov1 = cov_info1[ key ]
				except KeyError:
					cov1 = "-"
				try:
					cov2 = cov_info2[ key ]
				except KeyError:
					cov2 = "-"
				try:
					ho1, he1  = count_homo_and_hetero_vars( variants1[ key ] )
				except KeyError:
					ho1, he1  = 0, 0
				try:
					ho2, he2  = count_homo_and_hetero_vars( variants2[ key ] )
				except KeyError:
					ho2, he2  = 0, 0
				new_line = [ key, ho1, he1, cov1, ho2, he2, cov2, len( assembly[ key ] ) ]
				out.write( "\t".join( map( str, new_line ) ) + "\n" )
	
	# --- variant distribution plots --- #
	generate_var_distr_plots( variants1, variants2, assembly, min_contig_len, window_size, output_folder )
	
	# --- variant coverage analysis --- #
	try:
		analyze_data( frequencies1, frequencies2, coverages1, coverages2, output_folder )
	except ImportError:
		pass


if '--ref' in sys.argv and '--vcf1' in sys.argv and '--vcf2' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
