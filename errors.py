import sys
sys.stdout=open('Count of error types table.txt', 'w')#w means write, output will print to file instead of console, global 

import gzip

#LETTERS CODE to consider multiple letter changes in one error/diff
import re #using reg exp to find broken seq
letters=['A','C','T','G'] #creates a list
counter={'AC':0, 'CA':0,'AG':0, 'GA':0,'AT':0, 'TA':0,'CT':0, 'TC':0, 'CG':0, 'GC':0,'GT':0, 'TG':0} #creates a dictionary? of counters
ctable=[[],[],[]] #list of 3 lists/columns
ctable[0].append('Dif')
ctable[1].append('Con')
ctable[2].append('No.')
with gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt') as f:	
	for i, line in enumerate(f, 1): #contents is file as a string
	#TRYING TO RID INDEX OUT OF RANGE ERROR, try to get it to stop enumerating at end of file/last sequence
		if i % 4 == 2:
			seq=line #store the seq
			nextline= True #nextline, how specify, add 1 to i?
			if nextline: #inside the loop
				diffsline = gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt').readlines()[i] #reads file starting at 0
				diffs=diffsline.split() #makes a list from changes line, separated by spaces 
			nextline= False
			#extract pos no.
			for d in diffs: #and s in seqs to get them to align?
				pos=re.search(r'\d+',d) #search should get me first number found in each element, + sign means you have to find at least one
				if pos: #if you get something, to allow for no diffs
					#print (pos) #is a whole line with actual no. at end
					posno=int(pos.group())-1 #convert integer, manage offset
					cons=seq[posno] #that is cons #in the brackets is the indice of string
					for l in letters: #l is letter!
						if l != cons:                         
							if l in d: #look for l in a single diff
								counter[l+cons]+=1 
	for c in counter: #each counter is a line in the table
		#list of lists method
		splitc = re.findall('.', c) #splits c after a single position
		ctable[0].append(splitc[0]) #l is first component of c
		ctable[1].append(splitc[1]) #l and cons always G.. split c instead of retrieving after loop
		ctable[2].append(counter[c]) #list of ints? convert here to get list of strings 
	ctable[2]=[str(counter[c]) for counter[c] in ctable[2]] #or do this list comp
	#zip and joinrow method
	for row in zip(ctable[0],ctable[1],ctable[2]): #wants a str input #3rd element requires iteration? as many rows to match up with the first two arguments
		print('\t'.join(row)) #tab separated 

	#dataframe method to get headers, need to install pandas
	#import pandas
	#df=pandas.DataFrame(table[0],table[1],table[2])
	#df.columns["Diff","Cons","Count"]
	#print (df) 

#FLANK TABLE, sep file
sys.stdout=open('flank table.txt', 'w')
#loop through each seq and diff, go to error position in seq, return previous and after characters
letters=['A','C','T','G'] #creates a list
counter={'AC':0, 'CA':0,'AG':0, 'GA':0,'AT':0, 'TA':0,'CT':0, 'TC':0, 'CG':0, 'GC':0,'GT':0, 'TG':0} #creates a dictionary? of counters
ftable=[[],[],[],[],[]] #DIFFFFFFFFFFFFFFFFF
ftable[0].append('Seq')
ftable[1].append('Err')
ftable[2].append('Type')
ftable[3].append('Prev')
ftable[4].append('Post')
with gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt') as f:
	seqs=[] #initialise the list	
	for i, line in enumerate(f, 1):
		if i % 4 == 2:
			seq=line #store the seq
			print("seq:",seq)
			seqno=round((i/4)+0.1) #because 0.5 gets rounded down #DIFFFFFFFF
			print("seq no:",seqno)
			errorno={seqno:0} #right?
			nextline= True #nextline, how specify
			if nextline: #inside the loop
				diffsline = gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt').readlines()[i] #reads file starting at 0
				diffs=diffsline.split() #makes a list from changes line, separated by spaces 
				print ("diffs:",diffs)
			nextline= False
			for d in diffs: #and s in seqs to get them to align?
				pos=re.search(r'\d+',d) #search should get me first number found in each element, + sign means you have to find at least one
				if pos: #if you get something, to allow for no diffs
					posno=int(pos.group())-1 #convert integer, manage offset
					print("pos no:",posno)
					print(posno+2)
					print(len(seq))
					cons=seq[posno] #that is cons #in the brackets is the indice of string
					print("cons char:",cons)
					for l in letters: #l is letter!
						if l != cons:                         
							if l in d: #look for l in a single diff
								counter[l+cons]+=1 
								errorno[seqno]+=1 #counts diffs for each diffsline #DIFFFFFFFFF
								ftable[1].append(str(errorno))
								ftable[0].append(str(seqno))
								ftable[2].append(l+cons) #all lists of string
								#go back and get chars, shorter way??
								if posno==0: # == checks if two things have the same value
									prev2=' '
									post2=seq[posno+1]+seq[posno+2]
								elif posno==1:
									prev2=' '+seq[0]
									post2=seq[posno+1]+seq[posno+2]
								elif len(seq)-1<posno+1: #-1 to manage offset, make same as posno, start at 0
									post2=' '
									prev2=seq[posno-2]+seq[posno-1]	
								elif len(seq)-1<posno+2: 
									post2=seq[posno+1]+' '
									prev2=seq[posno-2]+seq[posno-1]
								else:
									prev2=seq[posno-2]+seq[posno-1]
									post2=seq[posno+1]+seq[posno+2]
								ftable[3].append(prev2)
								ftable[4].append(post2)
								print(ftable[0])
								print(ftable[1])
								print(ftable[2])
								print(ftable[3])
								print(ftable[4])
	#DIFFFFFFFFFFFF
	for row in zip(ftable[0],ftable[1],ftable[2],ftable[3],ftable[4]): #wants a str input #3rd element requires iteration? as many rows to match up with the first two arguments
		print('\t'.join(row)) #tab separated 
		#not printing anything??

#NEXT
#flanking
#run on whole file