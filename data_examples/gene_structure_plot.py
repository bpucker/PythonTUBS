### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### https://github.com/bpucker/APPLS ###
### https://www.cebitec.uni-bielefeld.de/~bpucker ###

__usage__ = """
			python gene_structure_plot.py\n
			--gff <GFF3_INPUT_FILE>\n
			--out <FULL_PATH_TO_OUTPUT_DIRECTORY>\n
			bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
			"""


import matplotlib.pyplot as plt
from operator import itemgetter
import re

# --- end of imports --- #

def construct_gene_structure_plot( gene_structures, outputfile ):
	"""! @brief constructs the plot of all elements of the gene structure """
	
	plt.close('all')
	fig = plt.figure( facecolor='white' )
	ax = fig.add_subplot( 111, frame_on=False )
	max_length = 0
	annotations = []
	for idx, gene_structure in enumerate( gene_structures ):
		
		if gene_structure['length'] > max_length:
			max_length = gene_structure['length']
		
		# --- set positions for plotting of the features of this gene --- #
		exon_base_lvl = idx*5+5
		intron_base_lvl = idx*5+4
		
		introns = []	#are constructed while parsing the exons
		
		annotations.append( gene_structure['annotation'] )
		
		
		# --- plotting all exons --- #
		for i, exon in enumerate( gene_structure['exons'] ):
			ax.plot( [ exon[0], exon[1] ], [ exon_base_lvl, exon_base_lvl ], linestyle='-', color='blue', linewidth=2.0 )
			if i < ( len( gene_structure['exons'] ) - 1 ):
				if not gene_structure['rv_strand']:	#if gene is on forward strand
					introns.append( ( exon[1], gene_structure['exons'][ i+1 ][0] ) )	#add new tuple to intron list
				else:	#if gene is on reverse strand
					introns.append( ( exon[0], gene_structure['exons'][ i+1 ][1] ) )	#add new tuple to intron list
		
		
		# --- plotting all the CDS parts --- #
		for CDS in gene_structure['CDS']:
			ax.plot( [ CDS[0], CDS[1] ], [ exon_base_lvl, exon_base_lvl ], linestyle='-', color='red', linewidth=4.0 )
		
		# --- plotting all introns --- #
		for intron in introns:
			mean = int( ( intron[0] + intron[1] ) / 2.0 )
			ax.plot( [ intron[0], mean, intron[1] ], [ exon_base_lvl, intron_base_lvl, exon_base_lvl ], linestyle='-', color='green', linewidth=1.0  )
	
	
	# --- add gene names --- #
	for idx, annotation in enumerate( annotations ):
		if gene_structures[idx]['rv_strand']:
			anno_col = 'blue'
		else:
			anno_col = 'green'
		ax.text( max_length, idx*5+5, "Col-0 gene"+gene_structure['annotation'], color=anno_col )
		
	
	# --- set axis and label xaxis --- #
	ax.set_ylim( 0, len( gene_structures )*6 + 6 )
	ax.set_xlim( 0, max_length+0.2*max_length )
	ax.set_xlabel( "position in gene [bp]" )
	
	# --- remove yaxis --- #
	ax.yaxis.set_visible( False )
	
	#plt.show()
	fig.savefig( outputfile, dpi=600 )
	plt.close('all')


def get_information_from_gff( gff_file, gene_IDs, pattern ):
	"""! @brief extacts information about one gene of interest from the gff file (all splice variants) """
	
	gene_structures = []
	
	with open( gff_file, "r" ) as f:
		line = f.readline()
		exons = []
		CDS = []
		length = 0
		prev_ID = ""
		while line:
			if line[0] != '#':
				try:
					gene_ID = re.findall( pattern, line )[0]
					if gene_ID in gene_IDs:
						if 'mRNA' in line or 'exon' in line or 'CDS' in line:
							if not prev_ID:
								prev_ID = gene_ID
							# --- save data of previous entry --- #
							if gene_ID != prev_ID:
								gene_structures.append( { 	'length': length,
															'exons': exons,
															'CDS': CDS,
															'annotation': prev_ID,
															'rv_strand': strand
														} )
								exons = []
								CDS = []
								prev_ID = gene_ID
								length = 0
							
							# --- collect data from line --- #
							parts = line.strip().split('\t')
							if parts[2] == 'exon':
								exons.append( ( int( parts[3] ), int( parts[4] ) ) )
							elif parts[2] == 'CDS':
								CDS.append( ( int( parts[3] ), int( parts[4] ) ) )
							elif parts[2] == 'mRNA':
								length = int( parts[4] ) - int( parts[3] )
							
							if '+' in parts[7]:
								strand = True
							else:
								strand = False
				except IndexError:
					pass
			line = f.readline()
		gene_structures.append( { 	'length': length,
									'exons': exons,
									'CDS': CDS,
									'annotation': prev_ID,
									'rv_strand': strand
								} )
	return gene_structures


def adjust_collected_gene_structure_values( gene_structures ):
	"""! @brief change setoff of gene structure to 0 via length of element """
	
	new_gene_structures = []
	
	print gene_structures
	for gene in gene_structures:
		all_exon_borders = []	#collect all exons borders positions to find maximal value
		for each in gene['exons']:	#iterating over all exons
			all_exon_borders.append( each[0] )	#and appending each values
			all_exon_borders.append( each[1] )
	
		if len( all_exon_borders ) > 0:
			length = max( all_exon_borders ) - min( all_exon_borders )
			diff = min( all_exon_borders )
			try:
				if all_exon_borders[2] < all_exon_borders[1]:
					rev_strand = True
				else:
					rev_strand = False
			except:
					rev_strand = False
			
			# --- calculate new exon positions --- #
			new_exons = []
			for exon in gene['exons']:
				new_exons.append( ( exon[0]-diff, exon[1]-diff ) )
			
			# --- calculate new CDS positions --- #
			new_CDS = []
			for part in gene['CDS']:
				new_CDS.append( ( part[0]-diff, part[1]-diff ) )
			
			# --- construct new gene structure --- #
			new_gene_structure = { 	'length': length,
									'exons': new_exons,
									'CDS': new_CDS,
									'annotation': gene['annotation'],
									'rv_strand': rev_strand		#True or False
								}
			# --- collect all new gene structures in new list --- #
			new_gene_structures.append( new_gene_structure )
		
	return new_gene_structures


def main( arguments ):
	
	araport_gff_file = arguments[ arguments.index( '--gff' )+1 ]
	prefix = arguments[ arguments.index( '--out' )+1 ]
	
	araport_pattern = "AT\dG\d{5}\.1"
	gene_IDs = [ "AT1G01000" ]
		
	
	# --- collect and adjust data for gene structure from GFF file of Araport11 --- #
	araport_gene_structures = get_information_from_gff( araport_gff_file, gene_IDs, araport_pattern )
	araport_gene_structures = adjust_collected_gene_structure_values( araport_gene_structures )
	
	
	# ---- drawing gene structure plots --- #
	for gene_structure in araport_gene_structures:
		output_file = prefix + gene_structure['annotation'] + ".png"
		construct_gene_structure_plot( [ gene_structure ], output_file )


if __name__ == '__main__':
	
	if '--gff' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
	print "all done!"
