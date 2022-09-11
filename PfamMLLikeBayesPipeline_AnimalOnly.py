import sys, re, datetime
from ete3 import Tree
import numpy as np
from numpy import random
from multiprocessing import Pool
from DatePfams import date_pfams
from LossLikelihood import likelihood_pipeline
from collections import defaultdict

### Using datasets generated using 'FollowNumbers.py': the animal distributions of only those Pfams found in one or more 
### animal species, which have appropriate z-scores.

### creates dictionaries of pfamIDs and species
def parse_pfam_file(pfamfile):
	pfam_dict = dict()
	for line in pfamfile.readlines():
		pfam_name = line.split(',')[0]
		pfam_dict[pfam_name] = set([re.sub('\n','',x) for x in line.split(',')[1:]])
	return pfam_dict
	
### loading in our initial pfam dictionary for the phylostratigraphy dataset (James et al.) and for interpro intergenic scan hits
with open("PfamUIDsTable_Animal_ex_contaminants.txt", 'r') as pfamfile, open("pfams_w_additions_Animal_ex_contaminants.txt", 'r') as interprofile:
	pfam_dict_original = parse_pfam_file(pfamfile)
	pfam_dict_interpro = parse_pfam_file(interprofile) 
	
### loading in our animal phylogenetic tree
tree_file_name = "animals.nwk"
tree = Tree(tree_file_name, format=1)


### initialisation dictionary: combining our interpro and cds pfam hits into a single dictionary 
### of lists, in order: total hits, interpro only hits. If pfam has no interpro hits, the list
### only contains the total hits.
updated_pfam_dict = defaultdict(list)

for pfam, species in pfam_dict_original.items(): 
	species = [x for x in species]
	interprospecies = pfam_dict_interpro[pfam]
	### species identified only from interpro hits:
	interprospecies = [x for x in interprospecies if x not in species]	
	### incorporate our 'interpro' species, the hits identified from scans of DNA previously identified as intergenic:
	comb_species = set(interprospecies+species)
	if len(comb_species) != 0:
		updated_pfam_dict[pfam] = (comb_species, interprospecies)
	### add in our unchanged pfams
	else:
		updated_pfam_dict[pfam] = set(species)


### initialisation parameters: false positive rate of original hits (p_false_pos_cds), and 
### interpro scan hits on intergenic data (p_false_pos_interpro)
parameters = [10**-6, 0.2]

### little wrapper function to pass multiple parameters to multiprocessing
def multi_arg_pool(partitioned_pfams):
	return likelihood_pipeline(partitioned_pfams, parameters, tree, updated_pfam_dict)

num_processes = 20
		
### chunking the pfam dataset into lists to run multiprocessing on
num_pfams = len(updated_pfam_dict)	
partitioned_pfams = []	
### The size of the chunk is an integer that's rounded from the number of pfams divided by the number of processes. Return a list of lists of chunked pfams
sizeOfChunk = int(num_pfams/num_processes)
partitioned_pfams = [list(updated_pfam_dict.keys())[i:i + sizeOfChunk] for i in range(0, num_pfams, sizeOfChunk)]


### some test lines to check the pipeline on a pfam subset			
#partitioned_pfams = [['PF17281'],['PF00014']]
#partitioned_pfams = ['PF00571','PF00572','PF00573','PF00574', 'PF10307']
#partitioned_pfams = [['PF00571','PF00572','PF00573','PF00574','PF00575','PF00576'], ['PF00578','PF00579','PF00581','PF00583','PF00584','PF00585']]

###this is how to run without the iterations, on a single pfam list, updating the likelihood (but also, no multiprocessing- which can be useful for quick testing)
# multi_arg_pool(partitioned_pfams)	

### starting the iterative process on all pfams, using a Bayesian procedure, recalculating the false positive rates and rerunning the ML pipeline
### using a maximum of 5 iterations

if __name__ == '__main__':

	for i in range(0, 5):

		res_list = []
		orig_list = []
	
		total_original = 0
		total_interpro = 0	
		excluded_original = 0
		excluded_interpro = 0
	
		new_pfam_dict = {}
	
		print('iteration: ' + str(i))
		print(' false positive parameters:')
		print(parameters)
# 		p = Pool(num_processes)
# 		ml_sp_search_out = p.map(multi_arg_pool, partitioned_pfams)
	

		if len(partitioned_pfams) == num_processes:
			p = Pool(num_processes)
			ml_sp_search_out = p.map(multi_arg_pool, partitioned_pfams)
		else:
			p = Pool(len(partitioned_pfams))
			ml_sp_search_out = p.map(multi_arg_pool, partitioned_pfams)
	
	
		for search_result in ml_sp_search_out:
			for pfam, results in search_result.items():
				res_list.append(results[1])
				orig_list.append(results[2])	
				new_pfam_dict[pfam] = results[0]	
			
	
		### compare new dictionary, new_pfam_dict with the initialisation dictionary, updated_pfam_dict (dictionary of lists, of order: total hits, interpro only hits)
	
		### recalculate false positive rates for original hits and interpro scan hits. Interpro hits are only counted as such
		### if they are not species identified by the original cds scans. 		
	
		for new_pfam, new_species in new_pfam_dict.items():
			
			interpro_only = updated_pfam_dict[new_pfam][1]
			original_only = [x for x in updated_pfam_dict[new_pfam][0] if x not in interpro_only]
					
			total_original = total_original + len(original_only)
			total_interpro = total_interpro + len(interpro_only)

			###	if the species composition has changed due to the false positive hit likelihood procedure		
			if len(updated_pfam_dict[new_pfam][0]) != len(new_species):
				updated_interpro = [x for x in new_species if x in interpro_only]
				updated_original = [x for x in new_species if x not in interpro_only]
			
				updated_original = len(original_only) - len(updated_original)
				updated_interpro = len(interpro_only) - len(updated_interpro)
				
				excluded_original = excluded_original + updated_original
				excluded_interpro = excluded_interpro + updated_interpro


		print('how many species do we exclude? From original: ' + str(excluded_original) + ' From interpro: '+str(excluded_interpro))	
		
		try:
			p_false_pos_cds = excluded_original/total_original
		except:
			p_false_pos_cds = 10**-6
		try:
			p_false_pos_interpro = excluded_interpro/total_interpro
		except:
			p_false_pos_interpro = 0.2
			
		print("p_false_pos_cds: " + str(p_false_pos_cds) + ", p_false_pos_interpro: " + str(p_false_pos_interpro) + '\n')				

		diff_in_parameters = abs(p_false_pos_cds - parameters[0]) + abs(p_false_pos_interpro - parameters[1])	
		
		parameters = [p_false_pos_cds, p_false_pos_interpro]

		with open('PfamMLLikeBayesPipeline_MLLoss_AnimalOnly_track.txt', 'a') as resfile, open('PfamMLLikeBayesPipeline_OrigLoss_AnimalOnly_track.txt', 'a') as origfile:
			resfile.write('iteration: ' + str(i) +', false positive parameters: ' + str(parameters[0])+', '+str(parameters[1])+'\n')

			resfile.write('\n'.join(res_list))
			resfile.write('\n')
			origfile.write('iteration: ' + str(i) + '\n')
			origfile.write('\n'.join(orig_list))            
			origfile.write('\n')	
	
		if diff_in_parameters < 10 ** -6:
			print('convergence reached! Final false positive probabilities:')
			print("p_false_pos_cds: " + str(p_false_pos_cds) + ", p_false_pos_interpro: " + str(p_false_pos_interpro) + '\n')
			break		
			
	
	with open('PfamMLLikeBayesPipeline_MLLoss_AnimalOnly.txt', 'wt') as resfile, open('PfamMLLikeBayesPipeline_OrigLoss_AnimalOnly.txt', 'wt') as origfile:
		resfile.write('\n'.join(res_list))
		origfile.write('\n'.join(orig_list))

