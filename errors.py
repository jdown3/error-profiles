import sys

filename=r'data\Jurkat_only_S5.consensus.100.fastq.gz' #what do about r, people won't input
#LETTERS CODE to consider multiple letter changes in one error/diff
letters=['A','C','T','G'] #creates a list
counter={'AC':0, 'CA':0,'AG':0, 'GA':0,'AT':0, 'TA':0,'CT':0, 'TC':0, 'CG':0, 'GC':0,'GT':0, 'TG':0} #creates a dictionary? of counters
ctable=[[],[],[]] #list of 3 lists/columns
ctable[0].append('Dif')
ctable[1].append('Con')
ctable[2].append('No.')

flankno=2  
ftable=[[],[],[],[],[],[],[]] 
ftable[0].append('Se#') 
ftable[1].append('Di#')
ftable[2].append('Er#')
ftable[3].append('Pre')
ftable[4].append('Dif')
ftable[5].append('Con')
ftable[6].append('Post')

from collections import defaultdict
from collections import OrderedDict
cfcounter=defaultdict(int)
cftable=[[],[],[],[],[]]
cftable[0].append('Pre')
cftable[1].append('Dif')
cftable[2].append('Con')
cftable[3].append('Pst')
cftable[4].append('No.')

import re #using reg exp to find broken seq
import gzip
with gzip.open(filename, 'rt') as f:	
	for i, line in enumerate(f, 1): #contents is file as a string
		if i % 4 == 2:
			seq=line.rstrip() #store the seq, gets rid of white space/line afterwards, \n is last char
			seqno=round((i/4)+0.1) #because 0.5 gets rounded down
			errorno={seqno:0}
			diffno={seqno:0} 
			nextline= True #nextline, how specify, add 1 to i?
			if nextline: #inside the loop
				diffsline = gzip.open(filename, 'rt').readlines()[i] #reads file starting at 0
				diffs=diffsline.split() #makes a list from changes line, separated by spaces 
			nextline= False
			for d in diffs: #and s in seqs to get them to align?
				diffno[seqno]+=1
				pos=re.search(r'\d+',d) #search should get me first number found in each element, + sign means you have to find at least one
				if pos: #if you get something, to allow for no diffs
					#print (pos) #is a whole line with actual no. at end
					posno=int(pos.group())-1 #convert integer, manage offset
					cons=seq[posno] #that is cons #in the brackets is the indice of string
					for l in letters: #l is letter!
						if l != cons:                         
							if l in d: #look for l in a single diff
								counter[l+cons]+=1 
								errorno[seqno]+=1 #counts diffs for each diffsline 
								ftable[2].append(str(errorno[seqno]))
								ftable[1].append(str(diffno[seqno]))
								ftable[0].append(str(seqno))
								ftable[4].append(l) #all lists of string
								ftable[5].append(cons)
								#go back and get chars, shorter way? slices
								if posno==0: # == checks if two things have the same value
									prev='  '
									post=seq[posno+1:(posno+1+flankno)] #slicing is inclusive of bounds, if flankno 2, like saying 1:2
								elif posno+1<=flankno: 
									prev=' '+seq[:posno]
									post=seq[posno+1:(posno+1+flankno)]
								elif len(seq)-1==posno: #-1 to manage offset, make same as posno, start at 0
									post='  ' #2 spaces 
									prev=seq[posno-flankno:posno]	
								elif posno+1>(len(seq)-flankno): #+1 to make it same as flankno, starting at 1
									post=seq[posno+1:]+' '
									prev=seq[posno-flankno:posno]
								else:
									prev=seq[posno-flankno:posno]
									post=seq[posno+1:(posno+1+flankno)]
								ftable[3].append(prev)
								ftable[6].append(post)
								with open('flank table.txt', 'w') as ft:
									for row in zip(ftable[0],ftable[1],ftable[2],ftable[3],ftable[4],ftable[5],ftable[6]): #wants a str input #3rd element requires iteration? as many rows to match up with the first two arguments
										ft.write('\t'.join(row)) #doesn't cumulative print
										ft.write('\n')
								#increment counter for that error seq
								cfcounter[prev+l+cons+post]+=1
								cfcounter[prev+l+cons+'  ']+=1
								cfcounter['  '+l+cons+post]+=1
								cfcounter[prev+' '+' '+post]+=1
								cfcounter[prev+' '+' '+'  ']+=1
								cfcounter['  '+' '+' '+post]+=1

with open('Count of error types table.txt', 'w') as ct: #w means write, output will print to file instead of console, global 
	for c in counter: #each counter is a line in the table
		#list of lists method
		splitc = re.findall('.', c) #splits c after a single position
		ctable[0].append(splitc[0]) #l is first component of c
		ctable[1].append(splitc[1]) #l and cons always G.. split c instead of retrieving after loop
		ctable[2].append(counter[c]) #list of ints? convert here to get list of strings 
	ctable[2]=[str(counter[c]) for counter[c] in ctable[2]] #or do this list comp
	#zip and joinrow method
	for row in zip(ctable[0],ctable[1],ctable[2]): #wants a str input #3rd element requires iteration? as many rows to match up with the first two arguments
		ct.write('\t'.join(row)) #tab separated 
		ct.write('\n')

	ct.write('\n')
	#count flanks table
	scfcounter=OrderedDict(reversed(sorted(cfcounter.items(), key=lambda t: t[1]))) #don't get, descending value
	for k in scfcounter: 
		splitk = re.findall('.',k) #len is 6 elements
		cftable[0].append(splitk[0]+splitk[1]) 
		cftable[1].append(splitk[2])
		cftable[2].append(splitk[3])
		cftable[3].append(splitk[4]+splitk[5])
		cftable[4].append(str(scfcounter[k]))
	for row in zip(cftable[0],cftable[1],cftable[2],cftable[3],cftable[4]):
		ct.write('\t'.join(row)) 
		ct.write('\n')

#NEXT
#change error and diff counters
#run on whole file
