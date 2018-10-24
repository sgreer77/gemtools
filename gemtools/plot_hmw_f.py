import os, sys, __main__ as main
import pandas as pd, numpy as np
import math
import rpy2.robjects as robj
import rpy2.robjects.pandas2ri # for dataframe conversion
from rpy2.robjects.packages import importr

def plot_hmw(outpre='out',**kwargs):

	if 'in_windows' in kwargs:
		infile= kwargs['in_windows']
	if 'out' in kwargs:
		outpre = kwargs['out']

	file_in="count_110.txt"
	out_plot="out.pdf"

	df=pd.read_table(infile,sep="\t")

	# First, create the barcode list
	bc_list = df.columns.tolist()[5:]

	# sort the barcodes according to where their reads map
	min_list=[]
	for bc in bc_list:
		df_tmp = df[['window_start',bc]]
		df_tmp = df_tmp.loc[df_tmp[bc]>0]
		min_val = df_tmp['window_start'].min()
		max_val = df_tmp['window_start'].max()
		min_list.append([bc,min_val,max_val])

	min_df = pd.DataFrame(min_list, columns=['bc','min_val','max_val'])
	min_df.sort_values(by='min_val', inplace=True)
	bc_order = min_df['bc'].tolist()

	bc_counter=1
	melt_list = []
	for bo in bc_order:
		df_bc = df[['window_start','window_end',bo]]
		df_bc = df_bc.loc[df_bc[bo]>0]
		if not df_bc.empty:
			df_bc['value']=bc_counter
			df_bc['variable']=bo
			df_bc=df_bc[['window_start','window_end','value','variable']]
			df_bc['bcs_split']= df_bc['variable'].apply(lambda x: x.split("-")[0])
			bc_counter=bc_counter+1
		melt_list.append(df_bc)

	m1 = pd.concat(melt_list)
	
	if m1.empty:
		print "No mappings to plot -- exiting -- start by checking that you are used the correct bam file to generate the input file"
		sys.exit()

	label_calc_low = m1['window_start'].min()
	label_calc_hi = m1['window_end'].max()

	label_y=math.ceil(m1['value'].max()/20)*20

	y_brk=50
	#m1_1<-subset(m1, window_start>label_calc_low & window_end<label_calc_hi)
	m1_1 = m1
	m1_1['window_start_div'] = m1_1['window_start']/1000000

	label_low=(round(label_calc_low/1000)*1000)/float(1000000)
	label_hi=(round(label_calc_hi/1000)*1000)/float(1000000)
	break_calc=(label_hi-label_low)/6

	# convert the pandas dataframe to an R dataframe
	robj.pandas2ri.activate()
	m1_r = robj.conversion.py2ri(m1_1)

	plotFunc = robj.r("""
	 library(ggplot2)
 
	function(df, label_low, label_hi, break_calc,label_calc_low,label_calc_hi, label_y, y_brk, outpre)	{
		theme_set(theme_bw(35))
		p1 <- ggplot(df,aes(x=(window_start/1000000),y=value)) +
			geom_point(size=0.6, colour="black") +
			xlab("chr10 (Mb)") +
			ylab("SV-specific barcode") +
			theme(plot.margin= unit(c(0.5, 1.5, 0.5, 0.1), "lines"),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				legend.position="none",
				panel.border = element_rect(colour = "black", fill=NA, size=2),
				axis.text=element_text(size=50),
				axis.text.x = element_text(size=20),
				axis.text.y = element_text(size=20),
				axis.title.y = element_text(size=25),
				axis.title.x=element_blank()) +
			scale_x_continuous(limits = c(label_calc_low/1000000, label_calc_hi/1000000), breaks=seq(label_low,label_hi,break_calc)) +
			scale_y_continuous(breaks=seq(0,label_y,y_brk))

		ggsave(outpre, plot=p1, width=15, height=4, units="in", dpi=600)
		}
	""")

	plotFunc(m1_r, label_low, label_hi, break_calc,label_calc_low,label_calc_hi, label_y, y_brk, outpre)


