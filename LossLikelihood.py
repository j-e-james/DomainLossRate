import sys, re, datetime
from ete3 import Tree
import numpy as np
from numpy import random
from multiprocessing import Pool
from DatePfams import date_pfams
from collections import defaultdict

verbose = False
verboseprint = print if verbose else lambda *a, **k: None

### returns the minium number of nodes from the basal clade which contain all species in the list sp, i.e the most likely 
### number of losses under parsimony
def get_monophyletic_nodes(species_list, basal_clade):
	node_species = [x.get_leaf_names() for x in basal_clade.iter_search_nodes()]
	### find the clades that contain only species that do not have the pfam.
	all_loss_nodes = [x for x in node_species if all(species in species_list for species in set(x))]
	all_loss_nodes.sort(key = len)
	### find the minimum number of clades that contain all species do not have the pfam
	loss_nodes = []
	flat_list = [item for sublist in loss_nodes for item in sublist]
	for loss_node in reversed(all_loss_nodes): 
		if all(sp not in flat_list for sp in loss_node):
			loss_nodes.append(loss_node)
			flat_list = [item for sublist in loss_nodes for item in sublist]
	### returns list of lists of species, with each node represented as descendent species.
	### the legth of the list is the minimum number of losses
	return(loss_nodes)	


### called iteratively by calc_ML_likeihood, calculates the maximum likelihood derivatives of the loss rate	
def calc_derivatives(loss_rate, sp_clade, losses, basal_clade, totalDistance):
	first_derivative = 0
	second_derivative = 0
	## use sp_clade tree, which contains only retained branches
	k = []
	for node in sp_clade.iter_descendants():
		if (float(1-np.exp(node.dist*loss_rate))) != 0:					
			k.append(float(node.dist) / float(1-np.exp(node.dist*loss_rate)))
	kprod = np.prod(k)
	
	log_like = -1 * loss_rate * totalDistance
	
	## use loss branches
	j = []				
	for clade in losses:
		distances = []
		loss_tree = basal_clade.copy()
		loss_tree.prune(clade)
		for node in loss_tree.iter_descendants():				
			distances.append(node.dist)		
		### of a clade which does not have the pfam, we only want to calculate loss likelihood 
		### on the most basal loss branch			
		if float(np.exp(max(distances)*loss_rate)) != 0:
			j.append(float(max(distances)) / float(np.exp(max(distances)*loss_rate)))

			first_derivative += (float(max(distances)) * np.exp(loss_rate * float(max(distances)))) / (np.exp(loss_rate * float(max(distances))) - 1)
			second_derivative += -np.exp(loss_rate * float(max(distances))) * (float(max(distances)) / (np.exp(loss_rate * float(max(distances))) - 1)) ** 2
			log_like += np.log(1-np.exp(-loss_rate*float(max(distances))))

	jprod = np.prod(j)
	initial_likelihood_lambda = (jprod*kprod)
	return initial_likelihood_lambda, log_like, first_derivative, second_derivative




def calc_ML_likelihood(pfam, tree, species):
	new_tree = tree.copy()
	basal_clade = new_tree.get_common_ancestor(species)		

	### if there are any losses:
	if len(basal_clade) != len(species):

		### calculate the minimum number of losses over the tree
		losses = get_monophyletic_nodes([x for x in basal_clade.get_leaf_names() if x not in species], basal_clade)
		print(pfam+ ' losses recorded')
		print(losses)
		### calculate the number of speciations over the tree (not counting speciations that occur
		### after a loss event)			
		node_species = [x.get_leaf_names() for x in basal_clade.iter_search_nodes()]
		node_species = [x for x in node_species if any(sp in species for sp in set(x))]			
		node_species.sort(key = len)
		nodes = []
		flat_list = []
		for node in reversed(node_species[1:]):
			flat_list = [item for sublist in nodes for item in sublist]
			if len(node) != 1 and len(node) != len(basal_clade):
				nodes.append(node)
			elif len(node) == 1:
				if node[0] not in flat_list:
					nodes.append(node)
		speciations = len(nodes)

		try:
			parsimony_loss = float(len(losses))/float(speciations)
		except:
			parsimony_loss = float(len(losses))
			print(pfam+' had a no speciation error')
		

		### now find loss rate, where time is calculated as the total time over the tree, prior to loss
		totalDistance = 0
		sp_clade = basal_clade.copy()
		### create tree containing only species with the pfam
		sp_clade.prune(species)
		for node in sp_clade.iter_descendants():
			totalDistance += node.dist	
	
		### we then calculate an initial loss rate for the pfam	
		loss_rate = float(len(losses))/totalDistance

		### initialise maximum likelihood method with the parsimony loss rate				
		new_loss_rate = loss_rate 
		tolerance = 10 ** (-7) ### 7 digit accuracy is desired
		epsilon = 10 ** (-14)  ### Don't want to divide by a number smaller than this
		max_iterations = 20    ### Don't allow the iterations to continue indefinitely
	
		### finding the zero of the derivative (first_derivative) using Newton's method
		for i in range(max_iterations):
			likelihood_lambda, log_like, first_derivative, second_derivative = calc_derivatives(new_loss_rate, sp_clade, losses, basal_clade, totalDistance)
			first_derivative = first_derivative - totalDistance
			old_loss_rate = new_loss_rate
			new_loss_rate = old_loss_rate - first_derivative/second_derivative

			if new_loss_rate > 1:
				likelihood_lambda, log_like, first_derivative, second_derivative = calc_derivatives(loss_rate, sp_clade, losses, basal_clade, totalDistance)
				new_loss_rate = loss_rate
				break
			elif abs(second_derivative) < epsilon:
				break
			### here is the updated version after talking with P. Nelson
			elif abs(first_derivative) <= tolerance:
				break
			elif np.isnan(new_loss_rate):
				likelihood_lambda, log_like, first_derivative, second_derivative = calc_derivatives(loss_rate, sp_clade, losses, basal_clade, totalDistance)
				new_loss_rate = loss_rate
				break
							
	else: 
		### NB- this will reset speciations counts to 0- but we don't need to count them, as this only 
		### kicks in if there are no losses
		losses = []
		likelihood_lambda = 0
		log_like = 0
		speciations = 0
		parsimony_loss = 0
		loss_rate = 0
		loss_rate_fraction = 0
		new_loss_rate = 0		
		totalDistance = 0
		basal_clade = tree.get_common_ancestor(species)
		for node in basal_clade.iter_descendants():
			totalDistance += node.dist	
		second_derivative = - 2 * totalDistance**2

	verboseprint(pfam+','+str(loss_rate)+','+str(new_loss_rate))
	
	date = date_pfams(species,tree)
	### throughout this code 'new_loss_rate' is the ML loss rate, while 'loss_rate' is the original parsimony loss rate.
	writestring = pfam+'\t'+str(len(species))+'\t'+str(speciations)+'\t'+str(len(losses))+'\t'+str(parsimony_loss)+'\t'+str(loss_rate)+'\t'+str(new_loss_rate)+'\t'+str(second_derivative)+'\t'+str(date)+'\t'+str(species)
	return(log_like, second_derivative, writestring)		






def likelihood_pipeline(pfam_list, parameters, tree, pipeline_start_pfam_dict):
	p_false_pos_cds = parameters[0]
	p_false_pos_interpro = parameters[1]
	
	updated_pfam_dict = {}
	
	for pfam in pfam_list:
	
		### run calc_ML_likelihood on original pfam species composition, including both original and interpro hits.
		### Mute the function after this function call to check for how this calcualted the initial estimate of loss rate and then maximises the likelihood.
		  
		original_likelihood, second_derivative, original_writestring = calc_ML_likelihood(pfam, tree, pipeline_start_pfam_dict[pfam][0])
		
		
		###now perform setup for false positive pfam hit removal: removing clades and assessing whether the likelihood improves.
		if original_likelihood != 0: ### if the original loss rate is not 0

			comb_species = list(pipeline_start_pfam_dict[pfam][0]) ### includes interpro and non-interpro hits
			interprospecies = list(pipeline_start_pfam_dict[pfam][1]) ### includes interpro hits only
		
			verboseprint(len(comb_species))
			best_likelihood = original_likelihood
			best_writestring = original_writestring
			species = comb_species		
			best_species_composition = comb_species

			ln_p_sp_list = []
		
			### identify largest number of possible monophyletic clades, i.e. 'A group composed of a collection of organisms, including the 
			### most recent common ancestor of all those organisms and all the descendants of that most recent common ancestor'
		
			tree_copy = tree.copy()	
			basal_clade = tree_copy.get_common_ancestor(comb_species)
			basal_clade.prune(comb_species)
		
			monophyletic_nodes = []
			for node_sp in basal_clade.iter_descendants():
				for node in tree.get_monophyletic(values=node_sp.get_leaf_names(), target_attr="name"):		
					monophyletic_nodes.append(node.get_leaf_names())

			### now we need identify the biggest node, and all other nodes required to get to our full complement of comb_species in the smallest possible 
			### number of monophyletic nodes
			monophyletic_nodes = sorted(monophyletic_nodes, key = len, reverse= True)
			sp_set = []	### checking when we have our full complement of species, so we can break out as soon as possible

			### largest monophyletic grouping
			largest = monophyletic_nodes[0]
			node_set = [largest]	
			for x in largest:
				sp_set.append(x)
		
			### now identify all other groups that add only new species
			for mono in monophyletic_nodes[1:]:
				if len(sp_set) == len(comb_species):
					break
				elif any(x in mono for x in sp_set):
					pass
				else:
					for x in mono:
						sp_set.append(x)
					node_set.append(mono)
		
		
			### we then exclude species in each node, after first counting the number of species excluded from the interpro and cds hits
			for node in node_set:
				cds_sp_excluded = [x for x in node if x not in interprospecies]
				interpro_sp_excluded = [x for x in node if x in interprospecies]
				newspecies = [x for x in comb_species if x not in node]
				### we then recalculate loss rates with the new species complement, 
				### unless the above steps have removed all species but 1
				if len(newspecies) > 1:
					verboseprint(pfam+' new species complement to test:' )
					verboseprint('species with pfam:'+ str(len(newspecies)))
				

					### find new likelihood, using different false positive likelihoods for cds and interpro hits
					ln_p_sp = len(newspecies) * np.log(1-p_false_pos_cds) + len(cds_sp_excluded) * np.log(p_false_pos_cds) + len(interpro_sp_excluded) * np.log(p_false_pos_interpro)

					if ln_p_sp in ln_p_sp_list:
						continue
					else:
						ln_p_sp_list.append(ln_p_sp)
								 
						log_like, second_derivative, writestring = calc_ML_likelihood(pfam, tree, newspecies)
				
						new_likelihood = log_like + ln_p_sp
	# 					print('new likelihood ' + str(new_likelihood))
	# 					print(len(newspecies))
	# 					print('best like so far ' + str(best_likelihood))

						if new_likelihood > best_likelihood or best_likelihood == 0:
							best_likelihood = new_likelihood
							best_species_composition = newspecies										
							best_writestring = writestring

			verboseprint('best likelihood: '+str(best_likelihood))
			verboseprint(len(best_species_composition))
			verboseprint(best_species_composition)
					
		else: ### if the loss rate is 0, we don't do the iterative flipping procedure, and set our species composition to the original dataset
			best_species_composition = list(pipeline_start_pfam_dict[pfam][0])	
			best_writestring = original_writestring			
					
					
		verboseprint(best_writestring)
		sys.stdout.flush()
		
		updated_pfam_dict[pfam] = best_species_composition, best_writestring, original_writestring

	return(updated_pfam_dict)




