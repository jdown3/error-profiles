import gzip
with gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt') as f: 
	contents=f.read() #contents is a string
	print (contents) 
	#opened gzip file, called it f, read and printed it to output

print ("number of times A occurs = %d" % (contents.count('A')))
print ("number of times C occurs = %d" % (contents.count('C')))
print ("number of times G occurs = %d" % (contents.count('G')))
print ("number of times T occurs = %d\n" % (contents.count('T'))) 
#number of times each nt occurs in whole file, not of use
#%d means the argument after the next pink % goes there and is a digit, doing this to print a mix of string and digit

with gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt') as f: #why do i have to open again, why/when does it close
	line=f.readlines()
	print (line[1])
	print ("number of times A occurs = %d" % (line[1].count('A')))
	print ("number of times C occurs = %d" % (line[1].count('C')))
	print ("number of times G occurs = %d" % (line[1].count('G')))
	print ("number of times T occurs = %d\n" % (line[1].count('T'))) 
	print (line[5])
	print ("number of times A occurs = %d" % (line[5].count('A')))
	print ("number of times C occurs = %d" % (line[5].count('C')))
	print ("number of times G occurs = %d" % (line[5].count('G')))
	print ("number of times T occurs = %d\n" % (line[5].count('T'))) 

#do this as a loop, where i is the line number, to count nts for every sequence in file
with gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt') as f: #why include
	for i, l in enumerate(f, 1): #i is index number for lines of file, l is each thing/line iterated in enumerate, calling the start of the file line 1
		if i % 4 == 2: #if index number get a remainder of 2 when divided be 4, so lines 2,6,10,14
			print ("Line number %d" % i)
			print ("Sequence number %d" % round(i/4)) #dvide index/line no by 4 to get sequence number, when python rounds 98.5 sometimes it rounds down
			print ("Sequence: \n %s" % l) #l is the actual line
			print ("Number of times A occurs = %d" % (l.count('A'))) 
			print ("Number of times C occurs = %d" % (l.count('C')))
			print ("Number of times G occurs = %d" % (l.count('G')))
			print ("Number of times T occurs = %d\n" % (l.count('T')))
#all working up to here

#now try counting combos
#count letters together in every sequence
			print ("Number of times AT occurs = %d" % (l.count('AT'))) #12 possible combos
			print ("Number of times CGC occurs = %d\n" % (l.count('CGC'))) 

#want to go through subsets of + line and see if two letters are there, 6 combos 
#TRYING TO GET TO 3RD LINE, NOT WORKING
#with gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt') as f: #why include
	#for i, l in enumerate(f, 1): #i is index number for lines of file, l is each thing/line iterated in enumerate, calling the start of the file line 1
		#if i>2:
			#print ("Line number %d" % i)
			#sequence number? 
			#INSERT re code on broken patterns
			#i+=4
		#else:
			#break #WHY DOESNT PRINT

#using reg exp to find broken seq
#FIGURE OUT
import sys
sys.stdout=open('output.txt', 'w')#w means write, output will print to file instead of console, global 
import re
with gzip.open(r'data\Jurkat_only_S5.consensus.100.fastq.gz', 'rt') as f: #why include
	for i, line in enumerate(f, 0):
		if i % 4 == 2: #starting at line 0 now
			print ("\nLine number %d" % i)
			print ("Sequence number %d" % round(i/4))
			counter={} #creates a dictionary? of counters
			#FIGURE OUT
			for match in re.finditer( r'A[CG\d]*C', line):#instead of just returning the first find, the for loop will iterate until it's gone through and found them all
				print ('\nAC found: ', match.group()) #how put new line at end??
				#increment counter and print at end?
				counter['AC']+=1
			for match in re.finditer( r'A[CG\d]*G', line):
				print ('\nAG found: ', match.group())
				counter['AG']+=1
			for match in re.finditer( r'A[CG\d]*T', line):
				print ('\nAT found: ', match.group()) 
				counter['AT']+=1
			#for match in re.finditer( r'C[CG\d]*A', line):
				#print ('\nCA found: ', match.group())
			for match in re.finditer( r'C[CG\d]*G', line):
				print ('\nCG found: ', match.group())
				counter['CG']+=1
			for match in re.finditer( r'C[CG\d]*T', line):
				print ('\nCT found: ', match.group())
				counter['CT']+=1	
			#for match in re.finditer( r'G[CG\d]*A', line):
				#print ('\nGA found: ', match.group())	
			#for match in re.finditer( r'G[CG\d]*C', line):
				#print ('\nGC found: ', match.group())	
			for match in re.finditer( r'G[CG\d]*T', line):
				print ('\nGT found: ', match.group())
				counter['GT']+=1	
			#for match in re.finditer( r'T[CG\d]*A', line):
				#print ('\nTA found: ', match.group())	
			#for match in re.finditer( r'T[CG\d]*C', line):
				#print ('\nTC found: ', match.group())
			#for match in re.finditer( r'T[CG\d]*G', line):
				#print ('\nTG found: ', match.group())
			for c in counter:
				print ("%s:%d" % (c, counter[c]))

#try counting from console not original file, idk how
#so write output to file and then count from that
#COUNTING TOTAL ERROR TYPES, less efficient method cause making a file
with open('output.txt') as o: #r means read
	changes=o.read()
	print ("\nThe number of A and C: %d" % (changes.count('AC')))
	print ("\nThe number of A and G: %d" % (changes.count('AG')))
	print ("\nThe number of A and T: %d" % (changes.count('AT')))
	print ("\nThe number of C and G: %d" % (changes.count('CG')))
	print ("\nThe number of C and T: %d" % (changes.count('CT')))
	print ("\nThe number of G and T: %d" % (changes.count('GT')))
	#what format for results? table?

#do better method, counter for each line, then sum (above)

#THURS
#allow for 4 letter changes
#flanking

#ADD Function
#def add(a,b): 
#	return a + b