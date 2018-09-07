import sys
import pandas as pd

def get_hmw_summary(**kwargs):
	
	if 'bc_counts' in kwargs:
		bc_c = kwargs['bc_counts']
	if 'out' in kwargs:
		outpre = kwargs['out']

	df = bc_c

	header = df.columns.tolist()
	bc_header = [b for b in header if b not in ["id","name","chrom","window_start","window_end"]]

	b_keep = []
	info_keep = []

	for b in bc_header:
		print b
		df_b = df[["id","name","chrom","window_start","window_end",b]]	
		df_b = df_b.loc[df[b]>0]
		if not df_b.empty:
			output = df_b.groupby(["id","name"], sort=False).agg({"window_start":min, "window_end":max, b: ['count','sum']}).reset_index()
			if len(output.index)<2:
				print "empty"
			else:
				print output
				output.columns = ['_'.join(col).strip() for col in output.columns.values]
				output.columns = ['id','name','end','start','num_windows','num_reads']
				output = output[['id','name','start','end','num_windows','num_reads']]
				#output['len'] = output['end'] - output['start']
				print output

		
				output = output.groupby('id', sort=False)['name','start','end','num_windows','num_reads'].apply(lambda x: pd.DataFrame(x.values)).unstack().reset_index()
				print output
				output.columns = ['id','name_1', 'name_2', 'start_1', 'start_2','end_1', 'end_2', 'num_windows_1', 'num_windows_2', 'num_reads_1', 'num_reads_2']
				output['len_1'] = output.apply(lambda row: int(row['end_1']) - int(row['start_1']), axis=1)
				output['len_2'] = output.apply(lambda row: int(row['end_2']) - int(row['start_2']), axis=1)
				output['mol_len'] = output.apply(lambda row: int(row['len_1']) + int(row['len_2']), axis=1)
				output['num_windows'] = output.apply(lambda row: int(row['num_windows_1']) + int(row['num_windows_2']), axis=1)
				output['num_reads'] = output.apply(lambda row: int(row['num_reads_1']) + int(row['num_reads_2']), axis=1)
				output['bc'] = b
		
				outlist = output.values.tolist()[0]
		
				print output
				print outlist
		

				b_keep.append(b)
				info_keep.append(outlist)	
		else:
			print "empty"

	keep_cols = ["id","name","chrom","window_start","window_end"] + b_keep
	final_df = df[keep_cols]
	final_df.to_csv("test" + ".bc_windows.subset.txt", sep="\t", index=False)


	info_df = pd.DataFrame(info_keep, columns=['id','name_1', 'name_2', 'start_1', 'start_2','end_1', 'end_2', 'num_windows_1', 'num_windows_2', 'num_reads_1', 'num_reads_2', 'len_1', 'len_2', 'mol_len', 'num_windows', 'num_reads', 'barcode'])

	info_df.to_csv(outpre, sep="\t", index=False)

	print len(final_df.columns)


