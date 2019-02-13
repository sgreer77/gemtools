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

	i=0
	for name, group in grouped:
		group_dedup = group.drop_duplicates(subset=['q_st','q_en'])
		
		if len(group_dedup.index)>1:
			i += 1
			assess_list.append([name,i])
		else:
			assess_list.append([name,"False"])
		
	print assess_list

	df_assess = pd.DataFrame(assess_list, columns=['read','interesting'])
	print df_assess

	df_comb = pd.merge(df, df_assess, on="read", how="left")

	#df_interesting = df.query('interesting == True')
	df_comb.to_csv(outfile, sep="\t", index=False)
