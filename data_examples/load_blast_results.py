### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### https://github.com/bpucker/APPLS ###
### https://www.cebitec.uni-bielefeld.de/~bpucker ###

def load_BLAST_results( input_file ):
	"""! @brief load all BLAST results from file """
	
	data = []
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			data.append( { 	'query': parts[0],
							'subject': parts[1],
							'query_start': int( parts[6] ),
							'query_end': int( parts[7] ),
							'score': float( parts[-1] ) 
						} )
			line = f.readline()
	return data
