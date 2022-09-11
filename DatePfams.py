import os, sys, json, csv, copy, datetime
from ete3 import Tree

'''
Assign age to a pfam 
General Age-Assigning Algorithm:
--------------------------------
The reason a different algorithm is used when only one species is present (algorithm [2]) is because the function [species_list].get_common_ancestor is used to find the subtree in the overall tree that has the common ancestor of all the species in the species list as the root node. If only one species is present, then the entire tree is returned as the subtree, invalidating the methodology used for algorithm [1] (half the length of the entire tree will always be assigned as the age of the Pfam in this case, so all Pfams with only one species associated with them will get the same, inaccurate age). 
-----
algorithm [1]: Used when more than one species is present:
   - The ete3 function [species_list].get_common_ancestor() is used to pull the subtree which has the common ancestor of all species in the species_list as the root node. 
   - The ete3 function subtree.iter_descendants("preorder") is used to iterate through the tree, starting at the root common node and going down one branch of the subtree. Because the Newick tree downloaded from TimeTree.org doesn't assign absolute ages to each node but rather relative ages, the length of each branch needs to be added to a total age variable. 
   - During the iteration process, it's checked to see whether the current position in the tree is a leaf or a node. If it's a node, we continue searching through the tree. As soon as a leaf is reached, the age of that leaf is added to the total age and we stop searching the tree. This is because we've completely gone down one branch of the subtree and have obtained the age. 
   - We then use the ete3 subtree.traverse("postorder") to find the age of the common ancestor in the tree. We divide this in half and add it to the total. This is the final age we upload to MySQL
algorithm [2]: Used when only one species is present
   - The algorithm for this case is much simpler than algorithm [1]. The species is located within the overall tree, the distance (age) is determined using the simple ete3 function node.dist and the age is divided in half. This is the age that is assigned to the Pfam
'''

# 
# 
# start_time = datetime.datetime.now()
# print('Script Executing\nCurrent Time: %s\n'%datetime.datetime.now())
# 
# # Data table where Pfams are stored
# PfamTable = 'loss_rates...'
# 
# # Newick tree used to date the pfams
# NewickTreeFilename = 'PhylogeneticTree_AllSpecies.nwk'
# 
# # First, the Newick tree is loaded
# t = Tree(NewickTreeFilename,format=1)
# 
# with open(PfamTable, 'r') as pfamResults:
# 	pfamResults = pfamResults.readlines()
# 	for pfam in pfamResults:
# 		speciesList = list(pfam.split('{')[1].split('}')[0])
# 
# 

### function takes a list of species, and the NewickTree	
def date_pfams(speciesList, TreeFile):
	# The list of species which have the pfam, after Bayesian curation
	numberOfSpecies = len(speciesList)

	# If more than one species exists for a particular pfam, then algorithm [1] is implemented (See: description in
	# the beginning of this script
	if numberOfSpecies != 1:
		# The subtree is pulled starting with the common ancestor of all species in the list as the root node
		subtree = TreeFile.get_common_ancestor(speciesList)
		# The total distance starts at zero and will be added to as each branch in the subtree is searched
		totalDistance = 0
		# We start at the base node and proceed down a branch of the subtree
		for node in subtree.iter_descendants("preorder"):
			# So long as our location is not a leaf, we add the relative age of the node to the total and continue to search the tree
			if node.is_leaf() == False:
				totalDistance += node.dist
			# Once we're located on a leaf, we have searched an entire branch of our subtree and, once we add the age of the leaf, we have acquired the
			# total age of our subtree (minus the age of the common ancestor)
			else:
				totalDistance += node.dist
				break
		# The variable maxLength is used to find the node that has the greatest number of children in the tree
		maxLength = 0
		# We then traverse the tree starting from the leaves and working our way toward the root node
		for node in subtree.traverse("postorder"):
			# We get the descendants of each node
			nodeDescendants = node.get_descendants()
			# We then compare the number of descendants to the greatest number found for a node thus far
			if len(nodeDescendants) > maxLength:
				# If it has more descendants than any thus far, we make a note of which node it was 
				rootNodeName = node.name
				# and find the age of that node
				rootNodeDistance = node.dist
		# The node with the greatest number of descendants in the common ancestor and we take its relative age and divide it in half (estimating the Pfam arose
		# halfway down that branch) and add it to the total
		totalDistance += (rootNodeDistance/2)
	else:
		# If only one species is associated with that Pfam, we implement algorithm [2]
		species =  TreeFile&speciesList[0]
		totalDistance = (species.dist)/2
	return totalDistance


