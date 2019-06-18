import mappy as mp

def align_contigs(**kwargs):

	if 'infile_fasta' in kwargs:
		infile = kwargs['infile_fasta']
	if 'out' in kwargs:
		outfile = kwargs['out']
	if 'genome' in kwargs:
		genome = kwargs['genome']
	if 'preset' in kwargs:
		preset = kwargs['preset']
	if 'nthreads' in kwargs:
		nthreads = kwargs['nthreads']

	a = mp.Aligner(str(genome), preset=preset, n_threads = nthreads)

	if not a: raise Exception("ERROR: failed to load/build index")

	outfile = open(outfile, 'w')

	outfile.write("read\tchr\tpos\tr_st\tr_en\tq_st\tq_en\tprimary\tstrand\tcs\tcigstr\tcigtup\n")

	for name, seq, qual in mp.fastx_read(infile):
		print name
		for hit in a.map(seq, cs=True):
			outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name, hit.ctg, hit.r_st, hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.is_primary, hit.strand, hit.cs, hit.cigar_str, hit.cigar))

	outfile.close()
