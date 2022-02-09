import sys
import os
from os import path
import re
import pandas as pd
import numpy as np
import time
import subprocess
import math as mth
import argparse

# argparse
if __name__ == '__main__':
	class CustomHelpFormatter(argparse.HelpFormatter):
		def _format_action_invocation(self, action):
			if not action.option_strings or action.nargs == 0:
				return super()._format_action_invocation(action)
			default = self._get_default_metavar_for_optional(action)
			args_string = self._format_args(action, default)
			return ', '.join(action.option_strings) + ' ' + args_string

	fmt = lambda prog: CustomHelpFormatter(prog)
	parser = argparse.ArgumentParser(description='Calculates the FineMAV scores', formatter_class=fmt)

	required = parser.add_argument_group('required flags')
	required.add_argument('-i', '--input-file', metavar='<file>', required=True, help='file containing the CHROM:POS, REF, ALT and AF (optional: AA, CADD_PHRED)')
	required.add_argument('-x', '--prefix', metavar='<prefix>', required=True, help='prefix assigned to the output files')
	required.add_argument('-r', '--reference-genome', metavar='<hg19|hg38>', required=True, help='hg19 or hg38')

	optional = parser.add_argument_group('optional flags')
	optional.add_argument('-v', '--vep-file', metavar='<file>', required=False, help='file containing the LOCATION with the AA and/or CADD_PHRED')
	optional.add_argument('-p','--penalty', metavar='<int|float>', required=False, help='penalty parameter')
	optional.add_argument('-c', '--chunksize', metavar='<int>', required=False, help='number of lines per chunk (default=200000)')
	args = parser.parse_args()

# determine the start time so that the "time taken" can be printed out in the log file
start_time = time.time()

input_file = args.input_file
vep_file = args.vep_file
prefix = args.prefix
reference_genome = args.reference_genome
input_penalty = args.penalty
chunksize = args.chunksize
temp_file = str('temporary_' + prefix + '.txt')

# get absolute path to resource (for PyInstaller)
def resource_path(relative_path):
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

# assign the data files (chrom.sizes) and wigToBigWig
# hg19_chrom was downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
hg19_chrom = resource_path("hg19.chrom.sizes")
# hg38_chrom was downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
hg38_chrom = resource_path("hg38.chrom.sizes")
# wigToBigWig was downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/wigToBigWig
wigToBigWig = resource_path("wigToBigWig")

# this header is used for finemav()
finemav_header = ''

# these are the compulsory columns needed to calculate the FineMAV scores (excluding AF)
compulsory_cols = ['LOCATION', 'ID', 'REF', 'ALT', 'AA', 'CADD_PHRED']

# check if input files exist
if path.exists(input_file) == False:
	print ('%s does not exist' % input_file)
	sys.exit(1)

# check if reference genome exists
if (reference_genome != "hg19") & (reference_genome != "hg38"):
	print ('Reference genome is not valid')
	sys.exit(1)

# determine the chunksize
c_size = 200000
if chunksize is not None:
	c_size = int(chunksize)

# check the input file first
input_header = pd.read_csv(input_file, nrows=0, sep='\t').columns.tolist()
input_header = [re.sub(r"#\s?", "", i) for i in input_header]
input_header = [re.sub(r"\[\d+\]", "", i) for i in input_header]

if input_header[0] == 'CHROM:POS':
	input_header[0] = 'LOCATION'

AF_list = [i for i in input_header if "AF_" in i]

popname_AF = []
for i in AF_list:
    i = re.sub("AF_", "", i)
    popname_AF.append(i)
popname_AF.sort()

print ("Populations detected: " + str(popname_AF)[1:-1])

# determine the penalty parameter
penalty_param = {2: 4.962,
                 3: 3.5,
                 4: 2.981,
                 5: 2.271,
                 6: 2.533,
                 7: 2.411}
penalty = 0
if input_penalty is None:
	if (len(popname_AF) > 7) | (len(popname_AF) < 2):
	    print ('The default penalty parameter is set if number of populations range from 2 to 7. Please use the --penalty flag to specify the penalty.')
	    sys.exit(1)
	else:
	    penalty = penalty_param[len(popname_AF)]
if input_penalty is not None:
	penalty = float(input_penalty)

### FUNCTIONS ###
# define function --- vep_check_and_merge
vepLine = ''
vep_header = ''

def vep_check_and_merge(input_file, vep_file):
	global vep_header
	global finemav_header
	global input_header
	global missing_cols

	stringToMatch = '#Location\t'
	with open (vep_file, 'r') as file:
		for line in file:
			if stringToMatch in line:
				vepLine = line

	vep_header = vepLine.lstrip('#').rstrip('\n')
	vep_header = vep_header.split("\t")
	for n, i in enumerate(vep_header):
		if i == 'Location':
			vep_header[n] = 'LOCATION'
		if i == 'Uploaded_variation':
			vep_header[n] = 'ID'

	# check if compulsory columns are there
	missing_cols = list(set(compulsory_cols) - set(input_header) - set(vep_header))
	if len(missing_cols) != 0:
		print('The following columns are missing from the VEP file: %s' % missing_cols)
		sys.exit(1)

	print('Merging the two files')
	for an_af_chunk in pd.read_csv(input_file, sep='\t', chunksize=c_size, header=0, names=input_header):
	    for aa_cadd_chunk in pd.read_csv(vep_file, chunksize=c_size, sep='\t', names=vep_header, comment='#'):
	        df = pd.merge(an_af_chunk, aa_cadd_chunk, how='inner', on=['LOCATION'])
	        df.to_csv('temporary_' + prefix + '.txt', mode='a', index=False, sep='\t', header=False)

	vep_header.remove("LOCATION")
	finemav_header = input_header + vep_header

# stats on the dataset
numOfSNPs = 0
numOfNoAncestral = 0
numOfAAisREF = 0
numOfAAisALT = 0
numOfAAisNeither = 0
numOfChunks = 0
currentChunk = 0

# create dictionary for maximum FineMAV scores
max_finemav_dict = dict.fromkeys(popname_AF)
for key in max_finemav_dict:
    max_finemav_dict[key] = 0

def stats(file):
	global numOfSNPs
	global numOfNoAncestral
	global numOfAAisREF
	global numOfAAisALT
	global numOfChunks
	global numOfAAisNeither

	for chunk in pd.read_csv(file, sep='\t', comment='#', chunksize=c_size, index_col=False, names=finemav_header, keep_default_na=False):
		numOfSNPs = numOfSNPs + (len(chunk))
		numOfNoAncestral = numOfNoAncestral + chunk['AA'].str.contains('N|\.|-', regex=True).sum()
		numOfAAisALT = numOfAAisALT + (chunk['AA'].str.lower() == chunk['ALT'].str.lower()).sum()
		numOfAAisREF = numOfAAisREF + (chunk['AA'].str.lower() == chunk['REF'].str.lower()).sum()
		numOfChunks = numOfChunks + 1

	numOfAAisNeither = numOfSNPs - numOfNoAncestral - numOfAAisREF - numOfAAisALT
	print('Total number of SNPs: %s' % numOfSNPs)


# define function -- finemav
def finemav(file):
	# figure out the FineMAV scores
	for chunk in pd.read_csv(file, sep='\t', comment='#', chunksize=c_size, index_col=False, names=finemav_header):

	    # determine the chunk number
	    global currentChunk
	    currentChunk = currentChunk + 1
	    print('Calculating the FineMAV scores (chunk %s of %s)' % (currentChunk, numOfChunks))

	    # rearrange columns
	    #global rearrange_header
	    chunk = chunk[compulsory_cols + AF_list]

	    # split LOCATION into CHROM and POS
	    chunk[['CHROM', 'POS']] = chunk['LOCATION'].str.split(":", expand=True)
	    chunk.drop(columns=['LOCATION'], inplace=True)
	    chunk['POS'] = chunk['POS'].astype(int)

	    # put CHROM and POS in the beginning
	    rearrange_header = chunk.columns.tolist()
	    rearrange_header = rearrange_header[-2:] + rearrange_header[:-2]
	    chunk = chunk[rearrange_header]

	    # ancestral allele
	    chunk["DER"] = chunk["AA"]

	    for pop in popname_AF:
	        # scenario 1
	        chunk["DER"] = chunk["ALT"]
	        chunk["DAF_" + pop] = chunk["AF_"+ pop]

	        # scenario 2
	        # need to set CADD to null if the derived is the reference allele
	        chunk.loc[chunk["AA"].str.lower() == chunk["ALT"].str.lower(), ["DER"]] = chunk["REF"]
	        chunk.loc[chunk["AA"].str.lower() == chunk["ALT"].str.lower(), ["DAF_" + pop]] = 1 - chunk["AF_" + pop]

	        # scenario 4
	        # chunk_AA_duplicate assigns DER = ALT
	        # this is a copy of the main chunk
	        chunk_AA_duplicate = chunk.loc[(chunk["AA"].str.lower() != chunk["REF"].str.lower()) & (chunk["AA"].str.lower() != chunk["ALT"].str.lower()) & (chunk["AA"] != 'N') & (chunk["AA"] != '.') & (chunk["AA"] != '-'),].copy()
	        chunk_AA_duplicate["DER"] = chunk_AA_duplicate["ALT"]

	    for pop in popname_AF:
	        chunk_AA_duplicate["DAF_" + pop] = chunk_AA_duplicate["AF_"+ pop]

	    for pop in popname_AF:
	        # scenario 4
	        # main chunk assigns DER = REF
	        chunk.loc[(chunk["AA"].str.lower() != chunk["REF"].str.lower()) & (chunk["AA"].str.lower() != chunk["ALT"].str.lower()), ["DER"]] = chunk["REF"]
	        chunk.loc[(chunk["AA"].str.lower() != chunk["REF"].str.lower()) & (chunk["AA"].str.lower() != chunk["ALT"].str.lower()), ["DAF_" + pop]] = 1 - chunk["AF_" + pop]

	        # scenario 3
	        chunk.loc[(chunk["AA"] == 'N') | (chunk["AA"] == '.') | (chunk["AA"] == '-'), ["DER"]] = "."
	        chunk.loc[(chunk["AA"] == 'N') | (chunk["AA"] == '.') | (chunk["AA"] == '-'), ["DAF_" + pop]] = 0

	    chunk = pd.concat([chunk, chunk_AA_duplicate], ignore_index=True)

	    chunk["SUM_OF_DER"] = float(0)

	    for pop in popname_AF:
	        chunk["SUM_OF_DER"] = chunk["SUM_OF_DER"] + chunk["DAF_" + pop]

	    chunk["DAP"] = float(0)

	    for pop in popname_AF:
	        chunk.loc[chunk["SUM_OF_DER"] != 0, 'DAP'] = chunk["DAP"] + ((chunk["DAF_" + pop]/chunk["SUM_OF_DER"])**penalty)

	    for pop in popname_AF:
	        chunk["FineMAV_" + pop] = chunk["CADD_PHRED"] * chunk["DAF_" + pop] * chunk["DAP"]
	        chunk.loc[chunk["CADD_PHRED"].isnull(), ["FineMAV_" + pop]] = chunk["DAF_" + pop] * chunk["DAP"]

	    # sort them in order
	    chunk.loc[(chunk["CHROM"] == "X"),["CHROM"]] = "23"
	    chunk.loc[(chunk["CHROM"] == "Y"),["CHROM"]] = "24"
	    chunk = chunk.astype({'CHROM': 'int'})
	    chunk.sort_values(['CHROM', 'POS'], inplace=True)
	    chunk = chunk.astype({'CHROM': 'object'})
	    chunk.loc[(chunk["CHROM"] == 23),["CHROM"]] = "X"
	    chunk.loc[(chunk["CHROM"] == 24),["CHROM"]] = "Y"
	    chunk.reset_index(drop=True, inplace=True)

	    # calculate maximum FineMAV score
	    for pop in popname_AF:
	        if max_finemav_dict[pop] < chunk["FineMAV_"+ pop].max():
	            max_finemav_dict[pop] = chunk["FineMAV_"+ pop].max()

	    # export table
	    print('Exporting the FineMAV scores into a table (chunk %s of %s)' % (currentChunk, numOfChunks))

	    if currentChunk == 1:
	        chunk.to_csv(prefix + '.txt', mode='a', sep='\t', index=False, na_rep=".", float_format='%.6f')
	    else:
	        chunk.to_csv(prefix + '.txt', header=False, mode='a', sep='\t', index=False, na_rep=".", float_format='%.6f')

	    # export to pre-bigWig
	    print('Exporting the FineMAV scores into temporary .wig files (chunk %s of %s)' % (currentChunk, numOfChunks))
	    for pop in popname_AF:
	        pop_max_finemav_df = chunk.loc[chunk.index.intersection(chunk.groupby(['CHROM', 'POS'])["FineMAV_" + pop].idxmax())]
	        chrom = pop_max_finemav_df["CHROM"].unique().tolist()

	        for i in chrom:
	            position = pop_max_finemav_df[pop_max_finemav_df["CHROM"] == i]["POS"]
	            scores = pop_max_finemav_df[pop_max_finemav_df["CHROM"] == i]["FineMAV_" + pop]
	            scores = pd.DataFrame({"#POS": position, "FineMAV_"+ pop: scores})
	            chrom_lines = '\nvariableStep chrom=chr' + str(i) + '\n\n{}'

	            # "with open" opens a specified file
	            # "a" is the appending mode which is used to add new data at the end of the file

	            with open('temporary_' + prefix + '_' + pop + '.wig', 'a') as f:
	                f.write(chrom_lines.format(scores.to_csv(index=False, sep='\t', float_format='%.6f')))

	# convert into bigWig
	print('Converting the .wig files into a bigWig files')

	if reference_genome == "hg19":
	    for pop in popname_AF:
	        subprocess.run([wigToBigWig, 'temporary_' + prefix + '_'+ pop + '.wig', hg19_chrom, prefix + '_' + pop + '.bw'])
	if reference_genome == "hg38":
	    for pop in popname_AF:
	        subprocess.run([wigToBigWig, 'temporary_' + prefix + '_' + pop + '.wig', hg38_chrom, prefix + '_' + pop + '.bw'])

	# delete the temporary files
	if vep_file is not None:
		os.remove('temporary_' + prefix + '.txt')
	for pop in popname_AF:
		os.remove('temporary_' + prefix + '_' + pop + '.wig')

	# creating the log file
	f = open(prefix +'_finemav.log', 'w')

	f.write('Name of input file: %s' % input_file)
	f.write('\nName of VEP file: %s' % vep_file)
	f.write('\nReference genome: %s' % reference_genome)
	f.write('\nChunksize: %s lines per chunk' % str(c_size))

	f.write('\n\nPenalty parameter: %s' % str(penalty))
	f.write('\nNumber of SNPs: %s' % str(numOfSNPs))

	f.write('\n\nNumber of SNPs where the ancestral allele...')
	f.write('\ndoes not exist: %s' % str(numOfNoAncestral))
	f.write('\nis the reference allele: %s' % str(numOfAAisREF))
	f.write('\nis the alternative allele: %s' % str(numOfAAisALT))
	f.write('\nis neither: %s' % str(numOfAAisNeither))

	f.write('\n\nNumber of populations: %s' % str(len(popname_AF)))
	f.write('\nPopulation names: %s' % ', '.join(popname_AF))

	for pop in popname_AF:
	    f.write('\nMaximum FineMAV score for %s: %s' % (pop, str('%.6f' % float(max_finemav_dict[pop]))))

	# write the time in float with 2 decimal places
	# create the time when the script ends
	end_time = time.time()
	f.write('\n\nTime taken: ' + str('%.2f' % float(end_time - start_time)) + ' seconds\nDone')

	f.close()

	print('Done')


# check if input_file contains compulsory_cols
missing_cols = list(set(compulsory_cols) - set(input_header))

if vep_file is None:
	if len(missing_cols) == 0:
	    finemav_header = input_header
	    stats(input_file)
	    finemav(input_file)
	else:
		print('You are missing the following columns: %s' % missing_cols)
		sys.exit(1)
else:
	vep_check_and_merge(input_file, vep_file)
	stats(temp_file)
	finemav(temp_file)
