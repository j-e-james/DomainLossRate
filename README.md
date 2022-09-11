# DomainLossRate
Scripts used in analysing rates of pfam domain loss rates across animal species.

pfams_w_additions_Animal_ex_contaminants.txt.gz
PfamUIDsTable_Animal_ex_contaminants.txt.gz
Pfam species hits, including only those Pfams identified in a minimum of 2 animal species in the dataset, and excluding Pfams with distributions indistinguishable from chance, using z-score method. 'PfamUIDsTable_Animal_ex_contaminants.txt.gz' are original Pfam hits, 'pfams_w_additions_Animal_ex_contaminants.txt.gz' include Pfam additional hits from scanning intergenic regions

animals.nwk: Phylogenetic tree of all animal species included in this analysis.

PfamMLLikeBayesPipeline_AnimalOnly.py: Main pipeline python script for analysing pfam loss rates.

LossLikelihood.py: Python modules called by main pipeline script.

DatePfams.py: Python modules called by main pipeline script.

Simulation.py: Python script used to generate simulated loss rates based on ISD and clustering values.

Rscript_Simulations.R: Analysis and plot generation of simulated data. Actual simulated values used are available as a table in this code.

FollowNumbers.py: Script written to perform a variety of follow the numbers tasks. Can also generate input files for further analysis steps.
Used to generate animal-specific datasets from the full Pfam dataset of James et al., and animal specific datasets for the intergenic genome scan data (aka interpro hits).  
Used to restrict dataset to only Pfams with z-scores < -2, indicating a random distribution across the total species phylogeny, either with or without the intergenic genome scan data (aka interpro hits).



