import sys
import argparse #for commandline interface, put in python, then directory to script
import re #using reg exp to find broken seq
import gzip
import contextlib #for exitstack
from collections import defaultdict
from operator import itemgetter #operator module has itemgetter for sorting output tables by value attribute

parser = argparse.ArgumentParser() #create an argument parser, need descript in brackets?\
parser.add_argument('file_direct', type=str, 
	help='what is the directory for the input file of sequences? Must be a .gz file.')
parser.add_argument('--flank', '-f', type=int, default=2, 
	help='how many nucleotides before and after the error position do you want to record? Default is 2.')
parser.add_argument('output_files_direct', type=str, 
	help='where do you want the output files to be stored and their name begin with?')
parser.add_argument('--output-flank', '-of', default=False, action='store_true', 
	help='whether you would like an output file of the flanking sequence of every error')
args = parser.parse_args()

filename = args.file_direct

#LETTERS CODE to consider multiple letter changes in one error/diff
letters = ['A','C','G','T'] #creates a list
#table to count each type of error, single letter change
counter = {'AC':0, 'CA':0,'AG':0, 'GA':0,'AT':0, 'TA':0,
         'CT':0, 'TC':0, 'CG':0, 'GC':0,'GT':0, 'TG':0} #creates a dictionary? of counters

acount = 0
ccount = 0
gcount = 0
tcount = 0

flankno = args.flank

diffs = []

cfcounter = defaultdict(int) #counting flanks
flcounter = defaultdict(int) #counting seg frags of flank length
seqsize = 2*flankno+1

with contextlib.ExitStack() as stack: #all files open in exit stack will be closed after with statement, future attempts to open won't work, keeps file open for length of with, won't be closed in next loop, so don't have to open again in every loop
	if args.output_flank:
		eft = stack.enter_context(open(args.output_files_direct + '.flankings' + '.txt', 'w'))
		eft.write('Seq Number\tDiff Number\tError Number\tPrev Seq\tDifference\tConsensus\tPost Seq\n')

	with gzip.open(filename, 'rt') as f:	
		for i, line in enumerate(f, 1):
			if i % 4 == 2:
				seq = line.rstrip() #store the seq, gets rid of white space/line afterwards, \n is last char
				acount+=(seq.count('A')) #working
				ccount+=(seq.count('C'))
				gcount+=(seq.count('G'))
				tcount+=(seq.count('T'))
				seqno = round(i/4 + 0.5) #because 0.5 gets rounded down
				print(seqno)
				startpos = 0
				while (len(seq)-startpos)>=seqsize: #as long as their are 5 chars left, length-startpos gets you how many chars left
					flcounter[seq[startpos:startpos+seqsize]]+=1
					startpos+=1
				errorno = 0
				diffno = 0
			if i % 4 == 3:	
				diffs = line.split() #makes a list from changes line, separated by spaces 
				for d in diffs: #after you get one diffsline, carry on
					diffno+=1
					pos = re.search(r'\d+',d) #search should get me first number found in each element, + sign means you have to find at least one
					if pos: #if you get something, to allow for no diffs
						posno = int(pos.group())-1 #convert integer, manage offset
						cons= seq[posno] #that is cons #in the brackets is the indice of string
						for l in letters: #l is letter!
							if l != cons:                         
								if l in d: #look for l in a single diff
									counter[l+cons]+=1 
									errorno+=1 #counts diffs for each diffsline 
									#go back and get chars, shorter way? slices
									if posno == 0: # == checks if two things have the same value
										prev = '  '
										post = seq[posno+1:(posno+1+flankno)] #slicing is inclusive of bounds, if flankno 2, like saying 1:2
									elif posno+1<=flankno: 
										prev = ' '+seq[:posno]
										post = seq[posno+1:(posno+1+flankno)]
									elif len(seq)-1 == posno: #-1 to manage offset, make same as posno, start at 0
										post = '  ' #2 spaces 
										prev = seq[posno-flankno:posno]	
									elif posno+1>(len(seq)-flankno): #+1 to make it same as flankno, starting at 1
										post = seq[posno+1:]+' '
										prev = seq[posno-flankno:posno]
									else:
										prev = seq[posno-flankno:posno]
										post = seq[posno+1:(posno+1+flankno)]
									if args.output_flank: #if you want a flanking table
											#write each current var. then forget
											eft.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seqno, diffno, errorno, prev, l, cons, post))
									#increment counter for that error seq
									cfcounter[prev+l+cons+post]+=1
									cfcounter[prev+l+cons+'  ']+=1
									cfcounter['  '+l+cons+post]+=1
									cfcounter[prev+' '+' '+post]+=1
									cfcounter[prev+' '+' '+'  ']+=1
									cfcounter['  '+' '+' '+post]+=1

with open(args.output_files_direct + '.other' + '.txt', 'w') as other:
	other.write('flankno\n%d' % (flankno))

with open(args.output_files_direct + '.nts' + '.txt', 'w') as nc: #w means write, output will print to file instead of console, global 
	nc.write('Nucleotide\tNumber\nA\t%d\nC\t%d\nG\t%d\nT\t%d\n' % (acount, ccount, gcount, tcount)) #dno how to order by value

with open(args.output_files_direct + '.errors' + '.txt', 'w') as ec: 
	ec.write('Difference\tConsensus\tNumber') #tab separated 
	scounter = sorted(counter.items(), key=itemgetter(1), reverse=True) #at this point only 2 cols in counter, e.g. AC 5
	#list of tuples, c is ('AA', 70)
	for c in scounter: #each counter is a line in the table
		splitc = re.findall('.', c[0]) #splits c after a single position
		ec.write('\n%s\t%s\t%d' % (splitc[0], splitc[1], c[1])) #list indice ,must be int or slice not str, but dict??

with open(args.output_files_direct + '.flankings_counts' + '.txt', 'w') as fc:
	#count flanks table
	fc.write('Prev Seq\tDifference\tConsensus\tPost Seq\tNumber')
	scfcounter = sorted(cfcounter.items(), key=itemgetter(1), reverse=True) #items sorts rows as one, reverse descending
	for k in scfcounter:  #len is 6 elements 
		fc.write('\n%s\t%s\t%s\t%s\t%d' % (k[0][0:flankno], k[0][flankno:flankno+1], 
			k[0][flankno+1:flankno+2], k[0][flankno+2:2*flankno+2], k[1]))
		#printing with extra notation

with open(args.output_files_direct + '.seq_fragments' + '.txt', 'w') as flc:
	flc.write('Fragment\tNumber')
	sflcounter = sorted(flcounter.items(), key=itemgetter(1), reverse=True)#just use sort, don't need ordereddict
	#sorts all items/rows together?
	for f in sflcounter: #f is tuple ('letters', 'no')
		flc.write('\n%s\t%d' % (f[0], f[1]))

#NEXT:
#analyse results in R, make standard program to analyse results, put together pipeline