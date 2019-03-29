import os
import sys
import __main__ as main
import pandas as pd
import numpy as np
import math
import rpy2.robjects as robj
import rpy2.robjects.pandas2ri # for dataframe conversion
from rpy2.robjects.packages import importr

def plot_hmw(outpre='out',**kwargs):

	if 'in_windows' in kwargs:
		infile= kwargs['in_windows']
	if 'out' in kwargs:
		outpre = kwargs['out']
	if 'sort_by_coord' in kwargs:
		sort_by_pos = kwargs['sort_by_coord']

	df=pd.read_table(infile,sep="\t")

	# remove any barcodes with no read mappings in the region
	bc_list_raw = df.columns.tolist()[5:]
	df_map = df.loc[(df[bc_list_raw].sum(axis=1) != 0), ]

	bc_list = df_map.columns.tolist()[5:]
	if len(bc_list)==0:
		print "No reads with these barcodes mapped in this region -- exiting"
		sys.exit()
	
	# count the number of chromosomes in the file -- can only plot 1 or 2
	chr_list = list(set(df_map['chrom'].tolist()))

	if len(chr_list)>2:
		print "Cannot plot breakpoints on more than 2 chromosomes at this time -- exiting"
		sys.exit()

	# if the user wants sorted barcodes, perform the sorting

	if sort_by_pos==False:
		bc_order = bc_list

	else: # if user wants sorting
		if len(chr_list)==1:
			chr_val = str(chr_list[0])

			min_list = []
			for bc in bc_list:
				df_tmp = df[['window_start',bc]]
				df_tmp = df_tmp.loc[df_tmp[bc]>0]
				if not df_tmp.empty:
					min_val = df_tmp['window_start'].min()
					max_val = df_tmp['window_start'].max()
					min_list.append([bc,min_val,max_val])

			min_df = pd.DataFrame(min_list, columns=['bc','min_val','max_val'])
			bc_order = min_df['bc'].tolist()
		
		elif len(chr_list)==2:

			global chr_val1
			global chr_val2
			chr_val1 = str(chr_list[0])
			chr_val2 = str(chr_list[1])


			df1 = df.loc[df['chrom']==chr_val1]
			df2 = df.loc[df['chrom']==chr_val2]

			min_list1=[]
			min_list2=[]        
			for bc in bc_list:
				df_tmp = df1[['window_start',bc]]
				df_tmp = df_tmp.loc[df_tmp[bc]>0]
				if not df_tmp.empty:
					min_val = df_tmp['window_start'].min()
					max_val = df_tmp['window_start'].max()
					min_list1.append([bc,min_val,max_val])
				else:
					df_tmp = df2[['window_start',bc]]
					df_tmp = df_tmp.loc[df_tmp[bc]>0]
					if not df_tmp.empty:
						min_val = df_tmp['window_start'].min()
						max_val = df_tmp['window_start'].max()
						min_list2.append([bc,min_val,max_val])

			min_df1 = pd.DataFrame(min_list1, columns=['bc','min_val','max_val'])
			min_df2 = pd.DataFrame(min_list2, columns=['bc','min_val','max_val'])

			if not min_df1.empty:
				min_df1.sort_values(by='min_val', inplace=True)
			if not min_df2.empty:
				min_df2.sort_values(by='min_val', inplace=True)
		
			min_df = pd.concat([min_df1,min_df2])
			bc_order = min_df['bc'].tolist()

	# melt the data frame in prep for plotting

	bc_counter=1
	melt_list = []
	for bo in bc_order:
		df_bc = df[['chrom','window_start','window_end',bo]]
		df_bc = df_bc.loc[df_bc[bo]>0]
		if not df_bc.empty:
			df_bc['value']=bc_counter
			df_bc['variable']=bo
			df_bc=df_bc[['chrom','window_start','window_end','value','variable']]
			df_bc['bcs_split']= df_bc['variable'].apply(lambda x: x.split("-")[0])
			bc_counter=bc_counter+1
			melt_list.append(df_bc)

	m1 = pd.concat(melt_list)


	if m1.empty:
		print "No mappings to plot -- exiting -- start by checking that you are used the correct bam file to generate the input file"
		sys.exit()
	
	plotFunc_1plot = robj.r("""
		suppressMessages(library(ggplot2)) 
		function(df,label_x,outpre)	{

		## general params
		#y_brk<-50
		label_y<-ceiling(max(df$value)/20)*20
		y_brk <- (ceiling((label_y/5)/10)*10)
	
		## calc df1/p1 params
		label_calc_low <- min(df$window_start)
		label_calc_hi <- max(df$window_end)
		label_low<-(round(label_calc_low/1000)*1000)/1000000
		label_hi<-(round(label_calc_hi/1000)*1000)/1000000
		break_calc<-(label_hi-label_low)/4
		df$window_start_div = df$window_start/1000000

		## plotting
		theme_set(theme_bw(35))
		p1 <- ggplot(df,aes(x=(window_start/1000000),y=value)) +
			geom_point(size=0.6, colour="black") +
			xlab(paste0(label_x," (Mb)")) +
			ylab("SV-spanning barcodes") +
			theme(plot.margin= unit(c(0.5, 1.5, 0.5, 0.1), "lines"),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				legend.position="none",
				panel.border = element_rect(colour = "black", fill=NA, size=2),
				axis.text=element_text(size=50),
				axis.text.x = element_text(size=20),
				axis.text.y = element_text(size=20),
				axis.title.y = element_text(size=25),
				axis.title.x=element_text(size=25)) +
			scale_x_continuous(limits = c(label_calc_low/1000000, label_calc_hi/1000000), breaks=seq(label_low,label_hi,break_calc)) +
			scale_y_continuous(breaks=seq(0,label_y,y_brk))
		#ggsave(outpre, plot=p1, width=15, height=4.5, units="in", dpi=100, device='png')
		ggsave(outpre, plot=p1, width=15, height=4.5, units="in", dpi=100)
		}
	""")

	plotFunc_2plot = robj.r("""
	 suppressMessages(library(ggplot2))
	 suppressMessages(library(gridExtra))
	function(df1, df2, label_x_1, label_x_2, outpre)	{
	
		## general params
		y_brk<-50
		df1_y<-ceiling(max(df1$value)/20)*20
		df2_y<-ceiling(max(df2$value)/20)*20
		label_y<-max(df1_y,df2_y)
		#y_brk <- (ceiling(max(label_y)/10)*10)/5
		y_brk <- (ceiling((label_y/5)/10)*10)
	
		## calc df1/p1 params
		label_calc_low_1 <- min(df1$window_start)
		label_calc_hi_1 <- max(df1$window_end)
		label_low_1<-(round(label_calc_low_1/1000)*1000)/1000000
		label_hi_1<-(round(label_calc_hi_1/1000)*1000)/1000000
		break_calc_1<-(label_hi_1-label_low_1)/4
		df1$window_start_div = df1$window_start/1000000

		## calc df2/p2 params
		label_calc_low_2 <- min(df2$window_start)
		label_calc_hi_2 <- max(df2$window_end)
		label_low_2<-(round(label_calc_low_2/1000)*1000)/1000000
		label_hi_2<-(round(label_calc_hi_2/1000)*1000)/1000000
		break_calc_2<-(label_hi_2-label_low_2)/4
		df2$window_start_div = df2$window_start/1000000

		## plotting
		theme_set(theme_bw(35))
		p1 <- ggplot(df1,aes(x=(window_start/1000000),y=value)) +
			geom_point(size=0.6, colour="black") +
			xlab(paste0(label_x_1," (Mb)")) +
			ylab("SV-spanning barcodes") +
			theme(plot.margin= unit(c(0.5, 1.5, 0.5, 0.1), "lines"),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				legend.position="none",
				panel.border = element_rect(colour = "black", fill=NA, size=2),
				axis.text=element_text(size=50),
				axis.text.x = element_text(size=20),
				axis.text.y = element_text(size=20),
				axis.title.y = element_text(size=25),
				axis.title.x=element_text(size=25)) +
			scale_x_continuous(limits = c(label_calc_low_1/1000000, label_calc_hi_1/1000000), breaks=seq(label_low_1,label_hi_1,break_calc_1)) +
			scale_y_continuous(breaks=seq(0,label_y,y_brk))
	
		p2 <- ggplot(df2,aes(x=(window_start/1000000),y=value)) +
			geom_point(size=0.6, colour="black") +
			xlab(paste0(label_x_2," (Mb)")) +
			ylab("SV-spanning barcode") +
			theme(plot.margin= unit(c(0.5, 1.5, 0.5, 2), "lines"),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				legend.position="none",
				panel.border = element_rect(colour = "black", fill=NA, size=2),
				axis.text=element_text(size=50),
				axis.text.x = element_text(size=20),
				axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				axis.title.x=element_text(size=25)) +
			scale_x_continuous(limits = c(label_calc_low_2/1000000, label_calc_hi_2/1000000), breaks=seq(label_low_2,label_hi_2,break_calc_2)) +
			scale_y_continuous(breaks=seq(0,label_y,y_brk))

		p_both<-grid.arrange(p1,p2,ncol=2,nrow=1, widths = c(1.1, 1))
	
		ggsave(outpre, plot=p_both, width=15, height=4.5, units="in", dpi=100)
		}
	""")

	# if there is one chromosome, plot like this:

	if len(chr_list)==1:
   
		robj.pandas2ri.activate()
		m1_r = robj.conversion.py2ri(m1)
	
		plotFunc_1plot(m1_r,chr_val,outpre)

	# if there are 2 chromosomes, plot like this:

	if len(chr_list)==2:
	
		m1_1 = m1.loc[m1['chrom']==chr_val1]
		m1_2 = m1.loc[m1['chrom']==chr_val2]
		
		robj.pandas2ri.activate()
		m1_1_r = robj.conversion.py2ri(m1_1)
		m1_2_r = robj.conversion.py2ri(m1_2)
	
		plotFunc_2plot(m1_1_r,m1_2_r,chr_val1,chr_val2,outpre)