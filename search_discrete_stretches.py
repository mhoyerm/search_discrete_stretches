import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Search for given peptide or nucleotide sequence in ORFs.', usage='[options]')
parser.add_argument('-F', '--fasta', type=str, required=True, help='fasta archive containing ORFs list')
parser.add_argument('--protein', default=False, action='store_true', help='use this command if the fasta archive is a protein file')
parser.add_argument('--keepframe', default=True, action='store_false', help='use this command if you want to respect the reading frame')
parser.add_argument('-E', '--peptides', type=str, default=False, help='list containing wanted peptides')
parser.add_argument('-N', '--nucleotides', type=str, default=False, help='list containing wanted nucletide sequence')
parser.add_argument('-O', '--output', type=str, required=True, help='output file (csv format)')

args = parser.parse_args()

def percent(currentgeneindex, genomelength): # Iterface updating the percentage of the ORFs read
	percentage = 100*currentgeneindex/genomelength
	interface = int(percentage)
	print str(interface) + '% complete\r',
import sys

def organizestretches(nucleotides, peptides):
	stretches = [[0,[],[]]]

	if nucleotides != False:
		nucleotides = nucleotides.replace('\r', '')
		nucleotides = str.split(nucleotides, '\n') # Split each line

		for currentcode in nucleotides:
			currentlength = len(currentcode)

			for stretchline in stretches:
				if stretchline[0] == 0:
					stretchline[0] = currentlength
					stretches.append([0,[],[]])
				if stretchline[0] == currentlength:
					stretchline[1].append(currentcode)
					break

	if peptides != False:
		peptides = peptides.replace('\r', '')
		peptides = str.split(peptides, '\n') # Split each line

		for currentcode in peptides:
			currentlength = len(currentcode)*3

			for i in range(len(stretches)):
				if stretches[i][0] == 0:
					stretches[i][0] = currentlength
					stretches.append([0,[],[]])
				if stretches[i][0] == currentlength:
					stretches[i][2].append(currentcode)
					break
	return stretches

def translate(orfsequence):
	aminoacids = {
		'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
		'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
		'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
		'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',

		'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
		'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
		'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
		'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',

		'TAT':'Y', 'TAC':'Y', 'TAA':'',  'TAG':'',
		'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
		'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
		'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',

		'TGT':'C', 'TGC':'C', 'TGA':'',  'TGG':'W',
		'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
		'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
		'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

	peptidesequence = ''
	for l in range(len(orfsequence)/3):
		codon = orfsequence[l*3] + orfsequence[l*3+1] + orfsequence[l*3+2]
		peptidesequence += aminoacids[codon]
	return peptidesequence

def evaluate(sequencecode, fout, geneid, stretcheslist):
	lencode = len(sequencecode)
	for j in range(lencode):
		if args.keepframe == True or args.protein == True or j % 3 == 0:
			for listline in stretcheslist:
				if listline[0] + j <= lencode:
					s = ''
					comparisonnucleotide = s.join(sequencecode[j : j+listline[0]])

					if args.nucleotides != False:
						for currentstretch in listline[1]:
							if currentstretch == comparisonnucleotide:
								relative = float(j) / float(lencode)
								fout.write("\n" + str(geneid) + ";" + str(j) + ";" + str(relative) + ";" + currentstretch)

					if args.peptides != False:
						if args.protein == False:
							comparisonpeptide = translate(sequencecode[j : j+listline[0]])
						if args.protein == True:
							comparisonpeptide = ''.join(sequencecode[j : j + listline[0]/3])
						for currentstretch in listline[2]:
							if currentstretch == comparisonpeptide:
								relative = float(j) / float(lencode)
								fout.write("\n" + str(geneid) + ";" + str(j) + ";" + str(relative) + ";" + currentstretch)
	return

def main():
	if args.protein == True and args.nucleotides != False:
		print 'Can not search for nucleotide sequence in a protein file'
		sys.exit()

	if args.protein == True and args.keepframe == True:
		print 'Can not ignore frame in a protein. Proceeding without ignoring it'

	file_in = open(args.fasta, 'r') # Input file
	in_file = file_in.read()

	try:
		file_peptides = open(args.peptides, 'r') # Input file
		peptides_file = file_peptides.read()
	except:
		peptides_file = ''

	try:
		file_nucleotides = open(args.nucleotides, 'r') # Input file
		nucleotides_file = file_nucleotides.read()
	except:
		nucleotides_file = ''

	out_file = open(args.output, 'w') # Output file
	
	actual_path = os.getcwd() # Get file path
	
	header_basic = 'Search for specific stretches of mRNA or polipeptide.\n'
	header_fasta = 'fasta file:;' + actual_path + '/' + args.fasta + '\n' # Prepare fasta path
	header_lists = 'peptide and nucleotide lists files:;' 
	if args.peptides != False:
		header_lists += str(args.peptides)
	if args.peptides != False and args.nucleotides != False:
		header_lists += ' | '
	if args.nucleotides != False:
		header_lists += str(args.nucleotides)
	header_lists += '\n'
	header_gene_stretch = '\n\n\n\n\n\n'
	header_data = 'name;position (in nucleotides);relative position;sequence'
	
	file_header = header_basic + header_fasta + header_lists + header_gene_stretch + header_data # Build header without difference value
	out_file.write(file_header) # Write header

	fastalist = in_file.replace('>UniRef100_','>')
	fastalist = str.split(fastalist, '>') # Split genes by '>'

	stretcheslist = organizestretches(nucleotides_file, peptides_file)

	for i in range(len(fastalist)): # Repeat for each gene
		percent(i, (len(fastalist)))

		# intrepret fasta
		sequencecode = ''
		fastalist[i] = fastalist[i].replace('\r', '')
		currentsequence = str.split(fastalist[i], '\n') # Split each line
		header = currentsequence[0] # The first line is the header
		geneid = str.split(header, ' ')[0]
		for j in range(1, len(currentsequence)):
			sequencecode += currentsequence[j] # Every line (except the first) is part of the code
		sequencecode = list(sequencecode) # List of all aminoacids in a gene

		evaluate(sequencecode, out_file, geneid, stretcheslist)
	print '\r100% complete!'

main()