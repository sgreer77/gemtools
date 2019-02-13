import pandas as pd

def assess_contigs(**kwargs):

	if 'infile_aln' in kwargs:
		infile = kwargs['infile_aln']
	if 'out' in kwargs:
		outfile = kwargs['out']

	df = pd.read_csv(infile, sep="\t", index_col = False)
	print df
	#df['interesting'] = False 
	#grouped = df.groupby('read').filter(lambda x: len(x)>=2) #only reads that appear 2+ times
	grouped = df.groupby('read')

	assess_list=[]

	j=0
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
		if i>1:
			if r_st_list[1] > r_en_list[0] or r_st_list[0] > r_en_list[1]: #interesting
				j += 1
				assess_list.append([name,j])
		else:
			assess_list.append([name,"False"])
		
	print assess_list

	df_assess = pd.DataFrame(assess_list, columns=['read','interesting'])
	print df_assess

	df_comb = pd.merge(df, df_assess, on="read", how="left")

	#df_interesting = df.query('interesting == True')
	df_comb.to_csv(outfile, sep="\t", index=False)
