import pandas as pd

def assess_contigs(**kwargs):

	if 'infile_aln' in kwargs:
		infile = kwargs['infile_aln']
	if 'out' in kwargs:
		outfile = kwargs['out']

	df = pd.read_table(infile, sep="\t", index_col = 0)
	df = df.reset_index() #otherwise 'read' column is index
	df['interesting'] = False 
	#grouped = df.groupby('read').filter(lambda x: len(x)>=2) #only reads that appear 2+ times
	grouped = df.groupby('read')

	for name, group in grouped:
		r_st_list = []
		r_en_list = []
		i = 0
		
		for row in group.iterrows():
			if abs(row[1]['r_en'] - row[1]['r_st']) > 500: # only consider mappings of at least 500
				r_st_list.append(row[1]['r_st'])
				r_en_list.append(row[1]['r_st'])
				i += 1
		print name
		print i
		if i == 2: #should generally have 2 rows
			if r_st_list[1] > r_en_list[0] or r_st_list[0] > r_en_list[1]: #interesting
				df.loc[name, 'interesting'] = True

	df_interesting = df.query('interesting == True')
	df_interesting.to_csv(outfile, sep="\t")
