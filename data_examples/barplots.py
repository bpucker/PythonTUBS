### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### https://github.com/bpucker/APPLS ###
### https://www.cebitec.uni-bielefeld.de/~bpucker ###

__usage__ = """
			python barplots.py\n
			--in <INPUT_VCF_FILE>\n
			--out <FULL_PATH_PREFIX_FOR_OUTPUT_FILES>\n
			feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
			"""


import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys, os

# --- end of imports --- #


def load_variants_from_vcf( vcf_file ):
	"""! @brief loads the variant informaiton from a VCF file """
	
	snps_per_chr = [ [], [], [], [], [] ]
	indels_per_chr = [ [], [], [], [], [] ]
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				try:
					parts = line.strip().split('\t')
					if len( parts[3] ) == len( parts[4] ):	#SNP
						snps_per_chr[ int( parts[0][-1] ) - 1 ].append( int( parts[1] ) )
					else:	#InDel
						indels_per_chr[ int( parts[0][-1] ) - 1 ].append( int( parts[1] ) )
				except:
					pass	#print line
			line = f.readline()
	return snps_per_chr, indels_per_chr


def construct_plot( snps_per_chr, indels_per_chr, chr_lengths, result_file, result_table, resolution, y_max_lim=4500 ):
	"""! @brief construct variant over Col-0 genome distribution plot """
	
	chromosome_length_units = [ 60, 40, 50, 40, 60 ]	#hard coded values!!
	
	f = plt.figure( figsize=(10,15) )
	gs = gridspec.GridSpec( 5, 1 )
	
	x_lim = [ 30, 20, 25, 20, 30 ]
	
	with open( result_table, "w" ) as out:
		for idx, chr_length in enumerate( chr_lengths ):
			chr_name = 'Chr' + str( idx + 1 )
			current_grid = gridspec.GridSpecFromSubplotSpec( 1, 60, subplot_spec = gs[ idx ] )
			
			# ---- identify positions of variants --- #
			upper_lim = resolution
			lower_lim = 0
			
			snp_data = []
			indel_data = []
			while True:
				if upper_lim >= chr_length:
					break
				else:
					snp_tmp = []
					indel_tmp = []
					for SNP in snps_per_chr[ idx ]:
						if SNP <= upper_lim and SNP > lower_lim:
							snp_tmp.append( 'X' )
					for indel in indels_per_chr[ idx ]:
						if indel <= upper_lim and indel > lower_lim:
							indel_tmp.append( 'X' )
					snp_data.append( len( snp_tmp ) )
					indel_data.append( len( indel_tmp ) )
				upper_lim += resolution
				lower_lim += resolution
			
			print "length of snp_data: " + str( len( snp_data ) )
			
			# --- improving x-axis --- #
			ax_a = plt.subplot( current_grid[ 0: chromosome_length_units[ idx ] ] )
			ax_a.set_xlim( 0, x_lim[ idx ] )
			ax_a.set_ylim( 0, y_max_lim )
			
			# --- plotting SNP and InDel distribution --- #
			for i, snps in enumerate( snp_data ):
				if i == 0:	#add labels only once!
					ax_a.plot( ( ((i*0.5)+0.2), ((i*0.5)+0.2) ), ( 0, snps ), "-", color="black", label="SNPs" )
					ax_a.plot( ( ((i*0.5)+0.3), ((i*0.5)+0.3) ), ( 0, indel_data[ i ] ), "-", color="red", label="InDels" )
				else:
					ax_a.plot( ( ((i*0.5)+0.2), ((i*0.5)+0.2) ), ( 0, snps ), "-", color="black" )
					ax_a.plot( ( ((i*0.5)+0.3), ((i*0.5)+0.3) ), ( 0, indel_data[ i ] ), "-", color="red" )
			
			# --- writing data into output table --- #
			out.write( 'Chr' + str( idx+1 ) + "SNPs:\t" + '\t'.join( map( str, snp_data ) ) + '\n' )
			out.write( 'Chr' + str( idx+1 ) + "InDels:\t" + '\t'.join( map( str, indel_data ) ) + '\n' )
			
			
			## --- improving ticks of y-axis a --- #
			max_yticks = 3
			yloc = plt.MaxNLocator(max_yticks)
			ax_a.yaxis.set_major_locator( yloc )
			
			
			current_labels = ax_a.get_xticks()
			labels = []
			for each in current_labels:
				labels.append( each  )
			ax_a.set_xticks([])
			ax_a.set_xticks( labels )
			
			ax_a.set_title( chr_name )
			ax_a.set_xlabel( "[ Mbp ]" )
			ax_a.set_ylabel( "number of variants" )	#"counts", "SNPs (black), InDel(red)"
			ax_a.legend( prop={'size':10} )
			
	gs.update( hspace=0.75 )
	plt.show()
	f.savefig( result_file, dpi=300 )


def main( arguments ):
	
	vcf_file = arguments[ arguments.index( '--in' )+1 ]
	prefix = arguments[ arguments.index( '--out' )+1 ]
	
	if prefix[-1] != '/':
		prefix += "/"
	
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	result_file = prefix + "genome_wide_variants.png"
	result_table = prefix + "genome_wide_variants.txt"
	
	resolution = 500000	#window size for analysis
	chr_lengths = [ 30427671, 19698289, 23459830, 18585056, 26975502 ]	#simplified script version !!
	snps_per_chr, indels_per_chr = load_variants_from_vcf( vcf_file )
	
	
	construct_plot( snps_per_chr, indels_per_chr, chr_lengths, result_file, result_table, resolution, y_max_lim=4500 )
	


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
	
	print "all done!"
