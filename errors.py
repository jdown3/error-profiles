import sys
sys.stdout=open('output.txt', 'w')#w means write, output will print to file instead of console, global 

import gzip
with gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt') as f: 
	contents=f.read() #contents is a string
	#print (contents) 
	#opened gzip file, called it f, read and printed it to output

#LETTERS CODE to consider multiple letter changes in one error/diff
import re #using reg exp to find broken seq
letters=['A','C','T','G'] #creates a list
counter={'AC':0, 'CA':0,'AG':0, 'GA':0,'AT':0, 'TA':0,'CT':0, 'TC':0, 'CG':0, 'GC':0,'GT':0, 'TG':0} #creates a dictionary? of counters
table=[[],[],[]] #list of 3 lists/columns
with gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt') as f:
	filestr=str(f)
	seqs=[] #initialise the list	
	for i, line in enumerate(f, 1): #contents is file as a string
	#TRYING TO RID INDEX OUT OF RANGE ERROR, try to get it to stop enumerating at end of file/last sequence
		if i % 4 == 2:
			seq=line #store the seq
			print ("\n\n\nSequence:",seq)
			nextline= True #nextline, how specify, add 1 to i?
			#continue #ignore logic in loop?
			if nextline: #inside the loop
				diffsline = gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt').readlines()[i] #reads file starting at 0
				print (diffsline)
				diffs=diffsline.split() #makes a list from changes line, separated by spaces 
				#do i want to store all the diff lines? 
				#print (diffs)
			nextline= False
			#extract pos no.
			for d in diffs: #and s in seqs to get them to align?
				###### got to allow for no diffs
				pos=re.search(r'\d+',d) #search should get me first number found in each element, + sign means you have to find at least one
				if pos: #if you get something
					#print (pos) #is a whole line with actual no. at end
					posno=int(pos.group())-1 #convert integer, manage offset
					print ('\nError location in consensus sequence: ', posno)
					print('\nConsensus character at error location: ', seq[posno])
					cons=seq[posno] #that is cons #in the brackets is the indice of string
					for l in letters: #l is letter!
						if l != cons:                         
							if l in d: #look for l in a single diff
								counter[l+cons]+=1 
			seqs.append(line) #adds seq line to list
			#print (seqs) #as a list
	for c in counter: #each counter is a line in the table
		print ("\nTotal %s errors:%d" % (c, counter[c])) #new total counts of error types, inclusive of multiple letter changes in one spot
		#list of lists(matrix) method
		splitc = re.findall('.', c)
		table[0].append(splitc[0]) #l is first component of c
		table[1].append(splitc[1]) #l and cons always G.. split c instead of retrieving after loop
		table[2].append(counter[c]) #list of ints?
	print(splitc[0])
	print(splitc[1])#need to be same format
	print(table[2])#don't need later
	table2=table[2]
	print (table2[0])
	#zip and joinrow method
	strtable2=str(table[2])
	print(strtable2)
	splittable2=strtable2.split(',')
	print(splittable2) #trying to get a list separated by square brackets, can't split
	print(splittable2[0])
	#ftable2=list(table2) 
	for row in zip(table[0],table[1],splittable2[0]): #wants a str input #3rd element requires iteration? do these arguments have a 1 char limit???
		print(' '.join(row)) #didn't print all, no borders

	#texttable method
	#from texttable import Texttable
	#t=Texttable()
	#t.addrows(['diff, 'cons', 'c', 'count'],[table])
	#print t.draw()
	#pink line & syntax error at cons

	#terminaltables method

	#asciitable method

	#tabulate method
	#from tabulate import tabulate
	#print tabulate ([table],[headers='diff','cons','c','count'],tablefmt='orgtb1')
	#got invalid syntax error after print tabulate 

	#one method
	#for line in table: #want to print the accumulation of errors after
		#print (' '.joinline())
			#print('\t'.joinline()) #\t is a tab 
			#error - has no joinline attribute

#zip method
#for row in zip(l + cons, counter[c]):
	#print (row)
			#error - zip arg #2 mjust support iteration


#NEXT
#fix sequence number, fix table/last col format , rid rest?
#flanking
