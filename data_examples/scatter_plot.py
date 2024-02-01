### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### https://github.com/bpucker/APPLS ###
### https://www.cebitec.uni-bielefeld.de/~bpucker ###

import matplotlib.pyplot as plt

fig, ax = plt.subplots( figsize=( 3, 3 ) )	#defining size of plot

x_values = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
y_values = [ 12, 6, 10, 4, 8, 14, 10, 11, 3, 9 ]

ax.scatter( x_values, y_values, color="red", s=50, marker="o", label="test" )
#setting color, marker size, marker shape and label of this group

ax.legend( numpoints=1 )
#each group is represented by only one marker in the legend (default=3)

ax.set_xlim( 0, 11 )	#set range of x-axis	
ax.set_ylim( 0, 15 )	#set range of y-axis

ax.set_xlabel( "pseudochromosome position [Mbp]" )

ax.spines["top"].set_visible(False)		#remove lines and ticks
ax.spines["right"].set_visible(False)	#remove lines and ticks

plt.subplots_adjust(left=0.1, right=0.99, top=0.97, bottom=0.2)
#adjust size of plot within figure

#plt.show()
fig.savefig( "my_plot.png", dpi=600 )	#write figure into output file
plt.close( "all" )	#destroy created figures (cleaning up)


