# DomainLossRate
Scripts used in analysing rates of pfam domain loss rates across animal species. Pipeline scripts were written by Jennifer James and Paul Nelson. Sara Willis wrote the original DataPfams script. All other python and R scripts written by Jennifer James.

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

Output.zip: Compressed folder of output files for PfamMLLikeBayesPipeline_AnimalOnly.py: 'PfamMLLikeBayesPipeline_MLLoss_AnimalOnly_track.txt', the ML loss rate results tracking changes through iterations, 'PfamMLLikeBayesPipeline_MLLoss_AnimalOnly.txt', which provides the final results ML loss rates and final species compositions as assigned through the false positive pipeline, 'PfamMLLikeBayesPipeline_OrigLoss_AnimalOnly_track.txt', parsimony loss rates tracking iterations (unchanging), 'PfamMLLikeBayesPipeline_OrigLoss_AnimalOnly.txt', the parsimony loss rates given the original species composition. 
Columns are: 'Pfam	Species	Speciations	Losses	ParsimonyLoss	ParsimonyLossRate	MLLossRate	SecondDerivative	Date	SpeciesList'
        



