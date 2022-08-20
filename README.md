# DomainLossRate
Scripts used in analysing rates of pfam domain loss rates across animal species.


Simulation.py: Python script used to generate simulated loss rates based on ISD and clustering values.

Rscript_Simulations.R: Analysis and plot generation of simulated data. Actual simulated values used are available as a table in this code.

FollowNumber.py: Script written to perform a variety of follow the numbers tasks. Can also generate input files for further analysis steps.
Used to generate animal-specific datasets from the full Pfam dataset of James et al., and animal specific datasets for the intergenic genome scan data (aka interpro hits).  
Used to restrict dataset to only Pfams with z-scores < -2, indicating a random distribution across the total species phylogeny, either with or without the intergenic genome scan data (aka interpro hits).

