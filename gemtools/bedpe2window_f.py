import sys
import pandas as pd

def get_types(y):
	y_list = y.split(";")
	y_dict = dict(s.split('=') for s in y_list)
	type = y_dict.get("TYPE")
	return type

def bedpe2window(**kwargs):
	if 'bedpe' in kwargs:
		sv_input = kwargs['bedpe']
	if 'window' in kwargs:
		wsize = kwargs['window']
	if 'out' in kwargs:
		outpre = kwargs['out']
	if 'mode' in kwargs:
		mode = kwargs['mode']
		
	mini_padder=200 #should be a little bigger than read length
	
	with open(sv_input) as f:
		for line in f:
			if line.startswith("#chr"):
				header_list = line.rstrip().split("\t")
			elif line.startswith('#'):
				pass
			else:
				break

	df = pd.read_table(sv_input, sep="\t", comment="#", header=None, names=header_list)
	
	coord_list = []

	if mode=="auto":
	
		df['type'] = df['info'].apply(lambda x: get_types(x))
	
		for index, row in df.iterrows():
			chrom1 = row['#chrom1']
			start1 = row['start1']
			stop1 = row['stop1']
			chrom2 = row['chrom2']
			start2 = row['start2']
			stop2 = row['stop2']
			sv_type = row['type']
			name = row['name']

			if sv_type=="DEL":
				wsize_del = wsize*2
				coord_list.append([chrom1,stop1-wsize_del,stop1,name,name,"in"])
				coord_list.append([chrom1,start2,start2+wsize_del,name,name,"in"])
				coord_list.append([chrom1,stop1+mini_padder,start2-mini_padder,name,name,"out"])

			elif sv_type=="DUP":
				wsize_check = wsize*2
				dup_size = start2 - stop1
				
				if dup_size<wsize_check:
					coord_list.append([chrom1,stop1,start2,name,name,"in"])
				else:
					coord_list.append([chrom1,stop1,stop1+wsize,name,name,"in"])
					coord_list.append([chrom1,start2-wsize,start2,name,name,"in"])	
						
				coord_list.append([chrom1,start1-wsize,start1-mini_padder,name,name,"out"])
				coord_list.append([chrom1,stop2+mini_padder,stop2+wsize,name,name,"out"])

			elif sv_type=="INV":
				sub_name_1 = str(name) + "_1"
				coord_list.append([chrom1,start1-wsize,start1,name,sub_name_1,"in"])
				coord_list.append([chrom1,start2-wsize,start2,name,sub_name_1,"in"])
				coord_list.append([chrom1,stop2+mini_padder,stop2+wsize,name,sub_name_1,"out"])

				sub_name_2 = str(name) + "_2"
				coord_list.append([chrom1,stop1,stop1+wsize,name,sub_name_2,"in"])
				coord_list.append([chrom1,stop2,stop2+wsize,name,sub_name_2,"in"])
				coord_list.append([chrom1,start1-wsize,start1-mini_padder,name,sub_name_2,"out"])		

			elif sv_type=="DISTAL":
				coord_list.append([chrom1,start1-wsize/2,stop1+wsize/2,name,name,"in"])
				coord_list.append([chrom2,start2-wsize/2,stop2+wsize/2,name,name,"in"])

			elif sv_type=="UNK":
				coord_list.append([chrom1,start1-wsize/2,stop1+wsize/2,name,name,"in"])
				coord_list.append([chrom2,start2-wsize/2,stop2+wsize/2,name,name,"in"])

			else:
				print "Unrecognized SV type in input"
				sys.exit()
				#out_list.append(list(row))
	
	elif mode=="window":
		
		for index, row in df.iterrows():
			chrom1 = row['#chrom1']
			start1 = row['start1']
			stop1 = row['stop1']
			chrom2 = row['chrom2']
			start2 = row['start2']
			stop2 = row['stop2']
			name = row['name']

			coord_list.append([chrom1,start1-wsize/2,stop1+wsize/2,name,name,"in"])
			coord_list.append([chrom2,start2-wsize/2,stop2+wsize/2,name,name,"in"])
	
	df_refined = pd.DataFrame(coord_list, columns = ['#chrom','start','stop','name','sub_name','status'])
	df_refined.to_csv(outpre, sep="\t", index=False)
	
	
	
