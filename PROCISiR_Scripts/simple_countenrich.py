import sys, os, math
import numpy as np
from scipy.interpolate import spline
from pylab import *

# aa_dict = {1:'A',2:'C',3:'D',4:'E',5:'F',6:'G',7:'H',8:'I',9:'K',10:'L',11:'M',12:'N',13:'P',14:'Q',15:'R',16:'S',17:'T',18:'V',19:'W',20:'Y',21:'*'}
#aa_dict = {'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19,'*':20}
aa_dict = {'K':0,'R':1,'E':2,'D':3,'Q':4,'N':5,'H':6,'W':7,'Y':8,'F':9,'M':10,'L':11,'I':12,'V':13,'A':14,'C':15,'S':16,'T':17,'P':18,'G':19,'*':20}

def all_indices(value, qlist):
    indices = []
    idx = -1
    while True:
        try:
            idx = qlist.index(value, idx+1)
            indices.append(idx)
        except ValueError:
            break
    return indices

def count_matrix(count_file,variable,ref_seq,pseudo):

	with open(count_file,'r') as initial_countfile:

		single_mutants = 0.0
		wt_counts = 0.0

		split_line = []

		for n,ini_line in enumerate(initial_countfile):
			split_line = ini_line.strip().split()

			if n==0:
				matrix = np.zeros([21,len(ref_seq)+1])
				matrix = matrix + pseudo

			elif n>0:
				if int(split_line[3]) == 1 and split_line[5] in aa_dict:
					matrix[aa_dict[split_line[5]],int(split_line[4])] += int(split_line[8])
					single_mutants += float(split_line[8]) + pseudo
				elif int(split_line[3]) == 0:
					wt_counts += float(split_line[8])
					single_mutants += float(split_line[8]) + pseudo

	variable_mask = (variable*np.ones([21,1])<1)
	masked_counts = np.ma.array(matrix,mask=variable_mask)
	shrunk_counts=np.ma.compress_cols(np.ma.masked_invalid(masked_counts))

	print single_mutants

	return(matrix,shrunk_counts,wt_counts/single_mutants)

def visualize(initial1,select1,variable1,ref_seq1,counts_cutoff,wt_enrich1,name):

	combined_initials = initial1
	combined_selected = select1

	ratios1 = np.subtract(np.log2(select1/np.sum(select1)),np.log2(initial1/np.sum(initial1)))
	ratios1 = np.subtract(ratios1,wt_enrich1)

	wt_enrich = (wt_enrich1)

	combined_ratios = ratios1

	n,bins,patches = hist(combined_ratios.flatten()[~np.isinf(combined_ratios.flatten())],bins=25,histtype='step')
	clf()

	bins = (bins[0:-1] + bins[1:])/2

	bnew = np.linspace(bins.min(),bins.max(),300)
	n_smooth = spline(bins,n,bnew)

	fig = figure()
	ax = fig.add_subplot(111)
	ax.scatter(combined_ratios,combined_initials,c=combined_ratios,cmap=cm.jet)
	ax.set_yscale('log')
	ylim(1,100000)

	ax.axhline(y=100,linestyle='--',c='k',lw=2)
	ax.axvline(x=0,linestyle='--',c='k',lw=2)

	line1 = ax.plot(bnew,n_smooth,c='k',lw=2)
	tick_params(top='off',right='off')
	xlabel('Log2 enrichment value',fontsize='large')
	ylabel('Log10 initial counts',fontsize='large')
	title(name+' Mutant enrichment vs. initial abundance',fontsize='large')

	savefig('enrichvsinit.png',dpi=500)
	clf()

	variable = variable1
	ref_seq = ref_seq1
	index_list = all_indices(1,variable)

	index_string = ''
	for x in index_list:
		index_string = index_string + str(x) + '\t'

	np.savetxt("counts.tsv",combined_initials,delimiter='\t',fmt='%1.1d')
	np.savetxt(name+"_sel_counts.tsv",combined_selected,delimiter='\t',fmt='%1.1d')
	np.savetxt(name+"_ratios.tsv",combined_ratios,delimiter='\t',fmt='%1.6f')#,header=index_string[0:-1])

	print sum(combined_initials)
	print sum(combined_selected)

	ratios_mask = (combined_initials < counts_cutoff)
	initials_mask = (np.zeros_like(combined_initials) > 1)

	for t in range(sum(variable)):
		initials_mask[aa_dict[ref_seq[index_list[t]]],t] = True
		ratios_mask[aa_dict[ref_seq[index_list[t]]],t] = True

	shrinkmasked_ratios = np.ma.array(combined_ratios,mask=ratios_mask)

	shrinkmasked_counts = np.ma.array(np.log10(combined_initials+1),mask=initials_mask)

	###Counts heatmap

        cmap=cm.winter

	imshow(shrinkmasked_counts,cmap=cmap,origin='lower',interpolation='none')
	tick_params(direction='out', top='off', right='off')
	tick_params(axis='y',labelleft='on')
	colorbar(shrink=.6)

	axvline(x=sum(variable)-.5,c='k',lw=1)

	for t in range(sum(variable)):
		annotate(ref_seq[index_list[t]],xy=(t-.4,aa_dict[ref_seq[index_list[t]]]-.4),xycoords='data',size='small')

	ylabel('AA at position')
	xlabel('Position in protein')
	# title('Log 10 counts in initial library')
	ylim(-.5,20.5)
	xlim(-0.5,sum(variable)-0.5)
#	yticks(np.arange(0,21,1),('K','R','E','D','Q','N','H','W','Y','F','M','L','I','V','A','C','S','T','P','G','*'),rotation=90,size='small')
        yticks(np.arange(0,21,1),('K','R','E','D','Q','N','H','W','Y','F','M','L','I','V','A','C','S','T','P','G','*'),rotation=90,size='small')
	xticks(np.arange(0,sum(variable)),[x+117 for x in all_indices(1,variable)],size=10,rotation=90)
	# show()

	fig = matplotlib.pyplot.gcf()
	fig.set_size_inches(20,5)

	savefig('combined_counts_plot.png',dpi=500)
	clf()

	###Ratios heatmap

        cmap=cm.coolwarm_r
        cmap.set_bad('0.75',1)

	imshow(shrinkmasked_ratios,cmap=cmap, origin='lower',interpolation='none',aspect='equal')
	tick_params(direction='out', top='off', right='off')
	tick_params(axis='y',labelleft='on')
	colorbar(shrink=.6)

	axvline(x=sum(variable1)-.5,c='k',lw=1)

	for t in range(sum(variable)):
		annotate(ref_seq[index_list[t]],xy=(t-.4,aa_dict[ref_seq[index_list[t]]]-.4),xycoords='data',size='small')

	ylabel('AA at position')
	xlabel('Position in protein')
	title(name+' Log 2 enrichment (initial counts >15)')
	xlim(-0.5,sum(variable)-0.5)
	ylim(-.5,20.5)
	# clim(-7.5,3.5)
	#yticks(np.arange(0,21,1),('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'),rotation=90,size='small')
	yticks(np.arange(0,21,1),('K','R','E','D','Q','N','H','W','Y','F','M','L','I','V','A','C','S','T','P','G','*'),rotation=90,size='small')
	xticks(np.arange(0,sum(variable)),[x+117 for x in all_indices(1,variable)],size=10,rotation=90)
	
	fig = matplotlib.pyplot.gcf()
	fig.set_size_inches(20,5)	

	savefig(name+'_combined_ratios_plot.png',dpi=500)
	# show()
	clf()

name = sys.argv[1]
ref_seq_1a = sys.argv[2]

variable_1a = ([1]*len(ref_seq_1a))
variable_1a.append(0)

baseline_counts1 = count_matrix(sys.argv[3],variable_1a,ref_seq_1a,0)
selected_counts1 = count_matrix(sys.argv[4],variable_1a,ref_seq_1a,1)

wt_enrich1 = np.log2((selected_counts1[2])/(baseline_counts1[2]))

visualize(baseline_counts1[1],selected_counts1[1],variable_1a,ref_seq_1a,15,wt_enrich1,name)
