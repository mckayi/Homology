import os
import sys
import re
import string, getopt, math, subprocess


####  Bcalm2Plex.py is a python script which takes in a *.fa dna file, creates a De Bruijn graph with
####  The kmers in the file, and outputs the vertices and edges of the graph in three text files.  
####  One text file lists the vertices and edges, another preps it for input in mathematica, and the third
####  Is a translator file, which lists the kmers and their id number for input in Plex or MMA. Ideally,
####  These .fa files have been run through Bcalm to obtain a "reduced" set of vertices for the graph.
####  Note that currently, only .fa files which started with 30mers are supported, but modifications to
####  This script will ensure that any size kmer will work.  The output .txt file is stored in a subfolder
####  Current_Directory\Output\inFileHandle, where inFileHandle is the .fa file without the extension.
####  For example, if G000464315.fna.bz2-30mers.bcalm.fa is your fasta file, then G000464315 is the folder
####  Where your text files are stored.  The .fa.txt file is ready to run through MATLAB, using Txt2ExplicitPlex.m
####  And the .fa.mma.txt is ready to run through Mathematica, in order to obtain higher dimensional coordinates.


	
def main(argv):
	filename = ''
	input_file = ''
	mathematica = 'Graph[{'
	mma = ''
	transfile = ''
	#Set up options
	try:
		opts, args = getopt.getopt(argv,"i:o:",["InputFile=","Outputfile="])
	except getopt.GetoptError:
		print('Unknown option, call using: ./LabelSortForPlex.py -i <InputFile.gexf> -o <Outputfile.txt>')
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i", "--InputFile"):
			input_file = arg
		elif opt in ("-o","--Outputfile"):
			filename = arg


	#Ensuring everything is saved to the correct path: Data is pulled from data folder, 
	#Output is stored in Output folder. If these folders don't exist, the script creates them.
	if filename == '':
		filename = input_file+'.txt'
		mma = input_file+'.mma.txt'
		tmpmma = input_file+'.tmp.txt'
		transfile = input_file+'.translator.txt'
	print("Input= " + input_file)
	print("Output= " + filename)
	path= os.path.abspath(__file__)
	path=path[0:len(path)-14]
	infileHandle = input_file[0:input_file.find('.fna')]
	outpath = path+'\\Output\\'+infileHandle
	datapath = path+'\\Data\\'+input_file
	if not os.path.exists(outpath):
		os.makedirs(outpath)

	filename = outpath+'\\'+filename
	mma = outpath+'\\'+mma
	tmpmma = outpath+'\\'+tmpmma
	transfile = outpath+'\\'+transfile
	#Remove Header
	my_re = re.compile(r">")
	f = open(datapath, 'r')
	print('Reading in vertices...')
	vertices = list()
	vertex_sizes = list()
	edges = list()
	for line in f:
		line_striped=line.strip()
		if not(my_re.search(line_striped)):
			vertices.append(line_striped)
		#else:
			#vertex_sizes.append(int(line_striped[1:]))
	f.close()
	print('Sorting Vertices...')
	#Sort Vertices, then store their position in the sorted list. Create two dictionaries:
	#Suffixes, and Prefixes, so that we may search and create our edges. The keys in the 
	#dicts are the prefixes (or suffixes), and the values are the index in the sorted list.
	#We also store the vertices in mathematica format so we may get higher dimensional coords.
	sortedVertices = sorted(set(vertices))
	prefixes = {}
	suffixes = {}
	f = open(filename, 'w')
	g = open(tmpmma,'w')
	translator = open(transfile,'w')
	g.write(mathematica)
	i=0
	for vertex in sortedVertices:
		i+=1
		translator.write(str(sortedVertices.index(vertex)) + '\t'+ vertex+'\n')
		if vertex[0:29] in prefixes:
			prefixes[vertex[0:29]].append(sortedVertices.index(vertex))
		else:
			prefixes[vertex[0:29]] = [sortedVertices.index(vertex)]
		if vertex[len(vertex)-30:] in suffixes:
			suffixes[vertex[len(vertex)-30:]].append(sortedVertices.index(vertex))
		else:
			suffixes[vertex[len(vertex)-30:]] = [sortedVertices.index(vertex)]
		f.write(str(sortedVertices.index(vertex)) + ', nan'+ '\n')
		if i<len(sortedVertices):
			g.write(str(sortedVertices.index(vertex))+', ')
		else:
			g.write(str(sortedVertices.index(vertex)))
	g.write('},{')
	translator.close()
	print('Forming Edges...')
	#Find edges, turn them into [x,y] format, where x and y are the positions of vertices in sorted list
	#Also write edges in Mathematica format as x->y
	k=0
	for prefix in prefixes:
		k+=1
		if "A"+prefix in suffixes:
			tempword = "A"+prefix
			for i in range(len(prefixes[prefix])):
				for j in range(len(suffixes[tempword])):
					edges.append(str(prefixes[prefix][i])+', '+str(suffixes[tempword][j]))
					g.write(str(prefixes[prefix][i])+'->'+str(suffixes[tempword][j])+', ')

		if "C"+prefix in suffixes:
			tempword = "C"+prefix
			for i in range(len(prefixes[prefix])):
				for j in range(len(suffixes[tempword])):
					edges.append(str(prefixes[prefix][i])+', '+str(suffixes[tempword][j]))
					g.write(str(prefixes[prefix][i])+'->'+str(suffixes[tempword][j])+', ')
		if "G"+prefix in suffixes:
			tempword = "G"+prefix
			for i in range(len(prefixes[prefix])):
				for j in range(len(suffixes[tempword])):
					edges.append(str(prefixes[prefix][i])+', '+str(suffixes[tempword][j]))
					g.write(str(prefixes[prefix][i])+'->'+str(suffixes[tempword][j])+', ')
		if "T"+prefix in suffixes:
			tempword = "T"+prefix
			for i in range(len(prefixes[prefix])):
				for j in range(len(suffixes[tempword])):
					edges.append(str(prefixes[prefix][i])+', '+str(suffixes[tempword][j]))
					g.write(str(prefixes[prefix][i])+'->'+str(suffixes[tempword][j])+', ')
	for edge in edges:
		f.write(edge + '\n')

	f.close()
	g.close()

	#Remove last comma and space from mathematica file, so that mathematica can read it in.
	g=open(tmpmma,'r')
	h=open(mma,'w')
	for line in g:
		h.write(line[0:len(line)-2])
	h.write("}, DirectedEdges->False]")
	g.close()
	h.close()
	#Clean up (remove temp mathematica file)
	os.remove(tmpmma)
	print("Complete!")


if __name__ == "__main__":
	main(sys.argv[1:])

