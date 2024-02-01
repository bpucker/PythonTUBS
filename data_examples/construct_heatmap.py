### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### https://github.com/bpucker/APPLS ###
### https://www.cebitec.uni-bielefeld.de/~bpucker ###


from flask import render_template, Flask
import re, pdfkit


app = Flask(__name__, template_folder="/PATH_TO_HTML_TEMPLATE/")


def load_gene_expression( gene_expression_file, seperator ):
	"""! @brief load column names, geneID and expression values from given file """
	
	# --- get max expression value --- #
	max_value = 0
	with open( gene_expression_file, "r" ) as f:
		header = f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split(seperator)
			tmp_max_val = max( map( float, parts[9:] ) )
			if tmp_max_val > max_value:
				max_value = tmp_max_val
			line = f.readline()
	
	
	# --- load gene expression values --- #
	gene_expression = []
	
	with open( gene_expression_file, "r" ) as f:
		column_names = [ "gene ID" ] + f.readline().strip().split(seperator)[9:]
		line = f.readline()
		while line:
			parts = line.strip().split(seperator)
			data = []
			for idx, part in enumerate( parts[9:] ):
				data.append( { 'title': parts[0] + '-' + column_names[ idx ], 'value': round( float( part ), 2), 'color': get_color( float( part ), max_value ) } )
			gene_expression.append( { 'name': parts[0], 'values': data } )
			line = f.readline()
	
	return column_names, gene_expression


def render_some_template( column_names, genes ):
	"""! @brief renders data into template and displays HTML document in browser via flask """
		
	return render_template( "heatmap.html", genes=genes, column_names=column_names )


def clamp( x ):
	"""! @brief checks for values outside of boundaries """
	return max( 0, min(x, 255) )


def get_color( value, max_value ):
	"""! @brief calculates background color based on current value and max_value """
	
	ratio = float( value ) / max_value
	r = 255 - int( 255. * ratio )			#define RGB color code
	g = 255		#define RGB color code
	b = 255						#define RGB color code
		
	hex_color = "#{0:02x}{1:02x}{2:02x}".format( clamp(r), clamp(g), clamp(b) )			#change RGB to hex color (for HTML)
	return hex_color



@app.route("/")
def main():
	"""! @brief used to call all other methods """
	
	# ---- file locations --- #
	gene_expression_data_file = "INPUT.csv"
	output_file = "heatmap.pdf"
	seperator=","
	
	# --- loading data --- #
	column_names, gene_expression = load_gene_expression( gene_expression_data_file, seperator )
	
	
	# --- construction of HTML template --- #
	rendered_template = render_some_template( column_names, gene_expression )
	
	
	# --- saving HTML to pdf --- #
	options = { 	'page-height': 400,	#CHANGE: dynamic size needed
					'page-width': 200
				}
	pdfkit.from_string( rendered_template, output_file, options=options )
	
	# --- returning the rendered HTML document for display via flask --- #
	return rendered_template


if __name__ == "__main__":
	app.debug = True	#leads to output of error messages on web page
	app.run()
	print "all done"
