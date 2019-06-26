#!/usr/bin/env python

import sys
import pandas as pd
import math
import rpy2.robjects as robj
import rpy2.robjects.pandas2ri # for dataframe conversion
from rpy2.robjects.packages import importr

def split_alleles(g):
	if "|" in g:
		a1 = g.split("|")[0]
		a2 = g.split("|")[1]
	elif "/" in g:
		a1 = g.split("/")[0]
		a2 = g.split("/")[1]
	else:
		a1 = "."
		a2 = "."
	return [a1,a2]


def plot_vars_and_blocks(**kwargs):

	if 'infile_basic' in kwargs:
		infile_basic = kwargs['infile_basic']
	if 'infile_blocks' in kwargs:
		infile_blocks = kwargs['infile_blocks']
	if 'region' in kwargs:
		region = kwargs['region']
	if 'out' in kwargs:
		plot_file_name = kwargs['out']


	df_basic = pd.read_table(infile_basic, sep="\t")
	df_blocks = pd.read_table(infile_blocks, sep="\t")
	chr = str(region.split(",")[0])
	start = int(region.split(",")[1])
	stop = int(region.split(",")[2])

	#start = int( math.floor(start_input / 100000.0) * 1000000.0 )
	#stop = int( math.ceil(stop_input / 1000000.0) * 1000000.0 )

	# get het SNVs in region (that pass filter)
	df_basic[['#chrom']] = df_basic[['#chrom']].astype(str)
	df_basic[['pos']] = df_basic[['pos']].astype(int)
	df_basic = df_basic[(df_basic['#chrom']==chr) & (df_basic['pos']>start) & (df_basic['pos']<stop) & (df_basic['var_type']=="snv") & (df_basic['filter']=="[]")]
	df_basic['a1'] = df_basic['gt'].apply(lambda x: split_alleles(x)[0])
	df_basic['a2'] = df_basic['gt'].apply(lambda x: split_alleles(x)[1])
	df_basic = df_basic.loc[df_basic['a1']!=df_basic['a2']]
	df_basic.rename(index=str, columns={"#chrom": "chrom"}, inplace=True)
	#df_basic.to_csv('basic.txt', sep="\t", index=False)

	# get phase blocks in region
	df_blocks[['chr']] = df_blocks[['chr']].astype(str)
	df_blocks[['beg_pos','end_pos']] = df_blocks[['beg_pos','end_pos']].astype(int)
	df_blocks = df_blocks.loc[(df_blocks['phased_het']>0) & (df_blocks['chr']==chr) & ( ((df_blocks['beg_pos']<start) & (df_blocks['end_pos']>stop)) | ((df_blocks['beg_pos']>start) & (df_blocks['beg_pos']<stop)) | ((df_blocks['end_pos']>start) & (df_blocks['end_pos']<stop)) )]
	df_blocks['beg_pos_check'] = df_blocks['beg_pos'].apply(lambda x: max(x,start))
	df_blocks['end_pos_check'] = df_blocks['end_pos'].apply(lambda x: min(x,stop)) 
	#df_blocks.to_csv('blocks.txt', sep="\t", index=False)

	# make variables for file name and x axis (to pass to R plotting function)
	#plot_file_name = str(infile_basic.split(".")[0]) + ".vars_blocks.png"

	if chr.startswith("chr"):
		x_axis_name = chr + " coordinate (Mb)"
	else:
		x_axis_name = "chr" + chr + " coordinate (Mb)"

	# now ready to plot
	plotFunc_blockvars = robj.r("""
		function(df_var,df_blk,b,e,xname,outplot)	{
			x_min = b/1000000
			x_max = e/1000000
			legend_x = x_min + (x_max-x_min)/2
			
			# adjust var table
			df_var$pos_mod = df_var$pos/1000000
			df_var$y = 1

			# adjust block table
			df_blk$beg_pos_plot<-df_blk$beg_pos_check/1000000
			df_blk$end_pos_plot<-df_blk$end_pos_check/1000000

			# open plot file
			png(outplot, width=9, height=3.2, units='in', res=300)
			par(mar=c(5,5,3,1)+.1)

			# make plot outline
			plot(y~pos_mod, data=df_var, type="n", xaxt="n", yaxt="n", ylab="", xlab="", xlim=c(x_min,x_max), ylim=c(0,3))
			#axis(1,at=x_lab,labels=x_lab, cex.axis=1.25,mgp=c(3,0.7,0))
			axis(1, cex.axis=1.25,mgp=c(3,0.7,0))
			mtext(side=1, cex=1.5, line=3, xname)
			axis(2,at=seq(1,2,1),labels=c("variants","blocks"), cex.axis=1.25, mgp=c(3,0.7,0), las=1)

			# plot variants
			points(y~pos_mod, data=df_var[df_var$phase_status=="not_phased",], pch=4, cex=1.6)
			points(y~pos_mod, data=df_var[df_var$phase_status=="phased",], pch=4, cex=1.6, col="red")

			# plot blocks
			rect(df_blk$beg_pos_plot, 1.8, df_blk$end_pos_plot, 2.2, lwd=2)

			# add legend
			par(xpd=TRUE)
			legend(legend_x, 4.5, c("phased", "unphased"), pch = c(4,4), col=c("red","black"), cex=1.4, pt.cex=1.6, horiz=T, bty="n", x.intersp=0.5)

			dev.off()
		}
	""")

	robj.pandas2ri.activate()
	df1 = robj.conversion.py2ri(df_basic)
	df2 = robj.conversion.py2ri(df_blocks)

	plotFunc_blockvars(df1,df2,start,stop,x_axis_name,plot_file_name)
