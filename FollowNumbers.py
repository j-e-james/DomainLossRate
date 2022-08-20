import sys, re, datetime
from ete3 import Tree
import numpy as np
from numpy import random
from multiprocessing import Pool
from DatePfams import date_pfams

######### Script written to perform a variety of follow the numbers tasks. Should not be run in one go.
######### Modules for specific tasks are separated by #########- the script was originally run in segments
######### by commenting out the rest- commenting out here does not indicate non-functional/unused code.

animal_list_master = set(["Ailuropoda_melanoleuca", "Anas_platyrhynchos", "Anolis_carolinensis", "Aotus_nancymaae", "Caenorhabditis_elegans", "Callithrix_jacchus", "Canis_lupus", "Capra_hircus", "Carlito_syrichta", "Cavia_aperea", "Cavia_porcellus", "Cebus_capucinus", "Cercocebus_atys", "Chinchilla_lanigera", "Chlorocebus_sabaeus", "Ciona_intestinalis", "Colobus_angolensis", "Cricetulus_griseus", "Cricetulus_griseus", "Danio_rerio", "Dasypus_novemcinctus", "Dipodomys_ordii", "Drosophila_melanogaster", "Eptatretus_burgeri", "Equus_caballus", "Felis_catus", "Ficedula_albicollis", "Fukomys_damarensis", "Gadus_morhua", "Gallus_gallus", "Gasterosteus_aculeatus", "Gorilla_gorilla", "Heterocephalus_glaber", "Homo_sapiens", "Ictidomys_tridecemlineatus", "Jaculus_jaculus", "Latimeria_chalumnae", "Lepisosteus_oculatus", "Loxodonta_africana", "Macaca_fascicularis", "Macaca_mulatta", "Macaca_nemestrina", "Mandrillus_leucophaeus", "Meleagris_gallopavo", "Mesocricetus_auratus", "Microcebus_murinus", "Microtus_ochrogaster", "Monodelphis_domestica", "Mus_caroli", "Mus_musculus", "Mus_pahari", "Mus_spretus", "Mustela_putorius", "Myotis_lucifugus", "Nannospalax_galili", "Nomascus_leucogenys", "Octodon_degus", "Ornithorhynchus_anatinus", "Oryctolagus_cuniculus", "Oryzias_latipes", "Otolemur_garnettii", "Ovis_aries", "Pan_paniscus", "Pan_troglodytes", "Panthera_pardus", "Panthera_tigris", "Papio_anubis", "Pelodiscus_sinensis", "Peromyscus_maniculatus", "Petromyzon_marinus", "Propithecus_coquereli", "Rattus_norvegicus", "Rhinopithecus_bieti", "Rhinopithecus_roxellana", "Saimiri_boliviensis", "Sarcophilus_harrisii", "Sus_scrofa", "Taeniopygia_guttata", "Takifugu_rubripes", "Tetraodon_nigroviridis", "Xenopus_tropicalis", "Xiphophorus_maculatus", "Acyrthosiphon_pisum", "Aedes_aegypti", "Amphimedon_queenslandica", "Anopheles_darlingi", "Anopheles_gambiae", "Apis_mellifera", "Atta_cephalotes", "Belgica_antarctica", "Bombyx_mori", "Brugia_malayi", "Caenorhabditis_brenneri", "Caenorhabditis_briggsae", "Caenorhabditis_japonica", "Caenorhabditis_remanei", "Capitella_teleta", "Crassostrea_gigas", "Culex_quinquefasciatus", "Danaus_plexippus", "Daphnia_pulex", "Dendroctonus_ponderosae", "Drosophila_ananassae", "Drosophila_erecta", "Drosophila_grimshawi", "Drosophila_mojavensis", "Drosophila_persimilis", "Drosophila_pseudoobscura", "Drosophila_sechellia", "Drosophila_simulans", "Drosophila_virilis", "Drosophila_willistoni", "Drosophila_yakuba", "Heliconius_melpomene", "Helobdella_robusta", "Lepeophtheirus_salmonis", "Lucilia_cuprina", "Mayetiola_destructor", "Melitaea_cinxia", "Mnemiopsis_leidyi", "Nematostella_vectensis", "Octopus_bimaculoides", "Pediculus_humanus", "Pristionchus_pacificus", "Rhodnius_prolixus", "Solenopsis_invicta", "Stegodyphus_mimosarum", "Strigamia_maritima", "Tribolium_castaneum", "Trichoplax_adhaerens", "Zootermopsis_nevadensis", "Acropora_digitifera", "Aedes_albopictus", "Apis_cerana", "Apis_dorsata", "Apis_florea", "Aplysia_californica", "Bactrocera_dorsalis", "Bactrocera_oleae", "Bemisia_tabaci", "Bicyclus_anynana", "Branchiostoma_belcheri", "Branchiostoma_floridae", "Camponotus_floridanus", "Ceratina_calcarata", "Ceratitis_capitata", "Ceratosolen_solmsi", "Cimex_lectularius", "Crassostrea_virginica", "Cyphomyrmex_costatus", "Drosophila_arizonae", "Drosophila_biarmipes", "Drosophila_bipectinata", "Drosophila_busckii", "Drosophila_elegans", "Drosophila_eugracilis", "Drosophila_ficusphila", "Drosophila_kikkawai", "Drosophila_miranda", "Drosophila_navojoa", "Drosophila_obscura", "Drosophila_rhopaloa", "Drosophila_serrata", "Drosophila_suzukii", "Drosophila_takahashii", "Dufourea_novaeangliae", "Echinococcus_granulosus", "Eufriesea_mexicana", "Eurytemora_affinis", "Metaseiulus_occidentalis", "Habropoda_laboriosa", "Hydra_vulgaris", "Limulus_polyphemus", "Linepithema_humile", "Mizuhopecten_yessoensis", "Monomorium_pharaonis", "Musca_domestica", "Myzus_persicae", "Nicrophorus_vespilloides", "Nilaparvata_lugens", "Onthophagus_taurus", "Montastraea_faveolata", "Orussus_abietinus", "Papilio_machaon", "Papilio_polytes", "Papilio_xuthus", "Parasteatoda_tepidariorum", "Pieris_rapae", "Pogonomyrmex_barbatus", "Priapulus_caudatus", "Pseudomyrmex_gracilis", "Saccoglossus_kowalevskii", "Stomoxys_calcitrans", "Trachymyrmex_cornetzi", "Trachymyrmex_septentrionalis", "Trachymyrmex_zeteki", "Varroa_destructor", "Vollenhovia_emeryi", "Wasmannia_auropunctata", "Bactrocera_cucurbitae", "Acinonyx_jubatus", "Balaenoptera_acutorostrata", "Bison_bison", "Bos_indicus", "Bos_mutus", "Bubalus_bubalis", "Camelus_bactrianus", "Camelus_dromedarius", "Camelus_ferus", "Castor_canadensis", "Ceratotherium_simum", "Chrysochloris_asiatica", "Condylura_cristata", "Echinops_telfairi", "Elephantulus_edwardii", "Eptesicus_fuscus", "Equus_asinus", "Equus_przewalskii", "Erinaceus_europaeus", "Hipposideros_armiger", "Leptonychotes_weddellii", "Lipotes_vexillifer", "Manis_javanica", "Marmota_marmota", "Meriones_unguiculatus", "Miniopterus_natalensis", "Myotis_brandtii", "Myotis_davidii", "Monachus_monachus", "Ochotona_princeps", "Odobenus_rosmarus", "Odocoileus_virginianus", "Orcinus_orca", "Orycteropus_afer", "Pantholops_hodgsonii", "Phascolarctos_cinereus", "Physeter_catodon", "Pongo_abelii", "Pteropus_alecto", "Pteropus_vampyrus", "Rhinolophus_sinicus", "Rousettus_aegyptiacus", "Sorex_araneus", "Trichechus_manatus", "Tursiops_truncatus", "Ursus_maritimus", "Vicugna_pacos", "Acanthisitta_chloris", "Acanthochromis_polyacanthus", "Alligator_mississippiensis", "Alligator_sinensis", "Anser_cygnoides", "Caprimulgus_carolinensis", "Apaloderma_vittatum", "Aptenodytes_forsteri", "Aquila_chrysaetos", "Balearica_regulorum", "Boleophthalmus_pectinirostris", "Buceros_rhinoceros", "Calidris_pugnax", "Callorhinchus_milii", "Calypte_anna", "Cariama_cristata", "Chaetura_pelagica", "Charadrius_vociferus", "Chelonia_mydas", "Chlamydotis_macqueenii", "Chrysemys_picta", "Clupea_harengus", "Colius_striatus", "Columba_livia", "Corvus_brachyrhynchos", "Coturnix_japonica", "Crocodylus_porosus", "Cuculus_canorus", "Cyprinodon_variegatus", "Cyprinus_carpio", "Egretta_garzetta", "Esox_lucius", "Eurypyga_helias", "Falco_cherrug", "Falco_peregrinus", "Fulmarus_glacialis", "Fundulus_heteroclitus", "Gavia_stellata", "Gavialis_gangeticus", "Gekko_japonicus", "Geospiza_fortis", "Haliaeetus_albicilla", "Haliaeetus_leucocephalus", "Haplochromis_burtoni", "Hippocampus_comes", "Ictalurus_punctatus", "Labrus_bergylta", "Larimichthys_crocea", "Lates_calcarifer", "Lepidothrix_coronata", "Leptosomus_discolor", "Lonchura_striata", "Manacus_vitellinus", "Maylandia_zebra", "Melopsittacus_undulatus", "Merops_nubicus", "Mesitornis_unicolor", "Nanorana_parkeri", "Neolamprologus_brichardi", "Nestor_notabilis", "Nipponia_nippon", "Nothobranchius_furzeri", "Notothenia_coriiceps", "Numida_meleagris", "Oncorhynchus_mykiss", "Opisthocomus_hoazin", "Parus_major", "Pelecanus_crispus", "Phaethon_lepturus", "Phalacrocorax_carbo", "Picoides_pubescens", "Poecilia_mexicana", "Pogona_vitticeps", "Protobothrops_mucrosquamatus", "Pseudopodoces_humilis", "Pterocles_gutturalis", "Pundamilia_nyererei", "Pygocentrus_nattereri", "Pygoscelis_adeliae", "Salmo_salar", "Scleropages_formosus", "Serinus_canaria", "Seriola_dumerili", "Sinocyclocheilus_anshuiensis", "Sinocyclocheilus_grahami", "Sinocyclocheilus_rhinocerous", "Stegastes_partitus", "Struthio_camelus", "Sturnus_vulgaris", "Tauraco_erythrolophus", "Thamnophis_sirtalis", "Tyto_alba", "Xenopus_laevis", "Zonotrichia_albicollis", "Ursus_arctos"])


######### as we are only calculating loss rates over animal species, we will restrict our input species to only 
######### animal species, and to only those pfams that are annotated in more than 1 animal.
# animalpfams = []
# with open("PfamUIDsTable_MV_Reformatted.txt", 'r') as pfamfile, open("PfamUIDsTable_Animal.txt", 'a') as animalfile:
# 	pfamfile = pfamfile.readlines()
# 	for line in pfamfile:
# 		pfam, species = line.split(',')[0], re.sub('\n','',line).split(',')[1:]
# 		new_species = [x for x in species if x in animal_list_master]
# 		if len(new_species) > 1:
# 			animalpfams.append(pfam)
# 			animalfile.write(pfam+','+','.join(new_species)+'\n') 
# 	
# with open("pfams_w_additions.txt", 'r') as pfamfile, open("pfams_w_additions_Animal.txt", 'a') as animalfile:
# 	pfamfile = pfamfile.readlines()
# 	for line in pfamfile:
# 		pfam, species = line.split(',')[0], re.sub('\n','',line).split(',')[1:]
# 		new_species = [x for x in species if x in animal_list_master]
# 		if len(new_species) >= 1 and pfam in animalpfams:
# 			animalfile.write(pfam+','+','.join(new_species)+'\n') 
	

######### compare interpro scan animal hits with the original animal hits
# with open("PfamUIDsTable_Animal.txt", 'r') as origfile, open("pfams_w_additions_Animal.txt", 'r') as interprofile, open("pfams_orig-interpro-differences_Animal.txt", 'a') as resfile:
# 	origfile = origfile.readlines()
# 	interprofile = interprofile.readlines()
# 	for line in interprofile:
# 		pfam, species = line.split(',')[0], re.sub('\n','',line).split(',')[1:]
# 		origpfam = [x for x in origfile if pfam in x]
# 		origspecies =  re.sub('\n','',origpfam[0]).split(',')[1:]
# 		diffs = list(set(species) - set(origspecies))
# 		if len(diffs) != 0:
# 			resfile.write(pfam+' '+str(len(diffs))+' '+str(len(origspecies))+' '+','.join(diffs)+'\n')
		
		
######### What Pfams should be excluded based on z-score values, from both original and interpro scan hits?
######### we exclude any Pfam that was both present in fewer than half the (animal) species in the dataset, and had a species distribution that was indistinguishable from chance
######### 'We exclude Pfams if they have z scores < -2 either with or without the intergenic hits flagged'- i.e. A and B
######### These values are taken directly from Paul Nelson- z-scores calculated using pfam_sigmas_add.txt
######### We exclude if Pfams turn up in either PfamsA (original) or PfamsB (interpro), so these pfams get removed from 
######### downstream analyses. This is not done iteratively- It's just one chunk of removing at the beginning, and so this
######### code is not rerun in the ML pipeline.
# with open("PfamsAExclusions.txt", 'r') as PfamsA, open("PfamsBExclusions.txt", 'r') as PfamsB, open("pfams_failing_zscore.txt", 'a') as resfile:	
# 	###121 in PfamsA, 220 in PfamsB, 80 Pfams in both= 261
# 	PfamsA= PfamsA.readlines()
# 	PfamsB = PfamsB.readlines()
# 	PfamsB = [re.sub('\n','',x) for x in PfamsB if x not in PfamsA]
# 	PfamsA = [re.sub('\n','',x) for x in PfamsA]
# 	PfamsFailingZscore = PfamsA + PfamsB 
# # 	for pfam in PfamsFailingZscore:
# # 		resfile.write(pfam) ### while this resultsfile is interesting, it is produced for all species, not only animals, and so length is misleading for this dataset
# 
######### Now we can generate our input file to the ML iterative process, which is restricted to animal species,
######### and excludes pfams with unacceptable z-scores.
### Load in the Pfamsfailingzscore list from above:
# with open("PfamUIDsTable_Animal.txt", 'r') as animalfile, open("pfams_w_additions_Animal.txt", 'r') as additfile, open("PfamUIDsTable_Animal_ex_contaminants.txt", 'a') as newanimalfile, open("pfams_w_additions_Animal_ex_contaminants.txt", 'a') as newadditfile:
# 	animalfile = animalfile.readlines()
# 	additfile = additfile.readlines()
# 	for line in animalfile:
# 		pfam = line.split(',')[0]
# 		if pfam not in PfamsFailingZscore:
# 			newanimalfile.write(line)
# 	for line in additfile:
# 		pfam = line.split(',')[0]
# 		if pfam not in PfamsFailingZscore:
# 			newadditfile.write(line)


with open("PfamUIDsTable_Animal_ex_contaminants.txt", 'r') as origfile:
	origfile = origfile.readlines()
	Pfams = [x.split(',')[0] for x in origfile]


with open('/Users/jennyjames/Desktop/LossRates/LossRatesForManuscript/PfamMLLikeBayesPipeline_MLLossAnimal_iteration4.txt', 'r') as LossFile, open('/Users/jennyjames/Desktop/LossRates/LossRatesForManuscript/PfamMLLikeBayesPipeline_MLLossAnimal_iteration4_ex_contaminants.txt', 'a') as resfile:
	LossFile = LossFile.readlines()
	for line in LossFile:
		pfam = line.split('\t')[0]
		if pfam in Pfams:
			resfile.write(line)	
 
	
	
	
	
	
### creates dictionaries of pfamIDs and species
def parse_pfam_file(pfamfile):
	pfam_dict = dict()
	for line in pfamfile.readlines():
		pfam_name = line.split(',')[0]
		pfam_dict[pfam_name] = set([re.sub('\n','',x) for x in line.split(',')[1:]])
	return pfam_dict	
### loading in our initial pfam dictionary for the phylostratigraphy dataset (James et al.) and for interpro intergenic scan hits
# with open("PfamUIDsTable_MV_Reformatted.txt", 'r') as pfamfile, open("pfams_w_additions.txt", 'r') as interprofile:
# 	pfam_dict_original = parse_pfam_file(pfamfile)
# 	pfam_dict_interpro = parse_pfam_file(interprofile)
