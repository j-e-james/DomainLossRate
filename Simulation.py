import numpy as np
from ete3 import Tree
"""

The non-monotonic relationship between loss rate and ISD remains statistically supported in the random effect model (p = 1e-7), with the model
 having a marginal R2, which represents the variance explained in loss rate just by the fixed effect of ISD, of 0.0051, and a slope of -2.42 
 for values of ISD below 0.18, and 1.19 for values above. The non-monotonic relationship with clustering also remained supported in the random 
 effect model (p = 0.017), with a marginal R2 0f 0.00056, with a slope of -0.29 for values of clustering below 0.85, and 0.099 for values of 
 clustering above 0.85.
 
"""


### loading in our overall phylogenetic tree, as used in James et al.
tree_file_name = "PhylogeneticTree_AllSpecies.nwk"
tree = Tree(tree_file_name, format=1)
### tree is 457 species

#tree.prune(["Ailuropoda_melanoleuca", "Anas_platyrhynchos", "Anolis_carolinensis", "Aotus_nancymaae", "Caenorhabditis_elegans", "Callithrix_jacchus", "Canis_lupus", "Capra_hircus", "Carlito_syrichta", "Cavia_aperea", "Cavia_porcellus", "Cebus_capucinus", "Cercocebus_atys", "Chinchilla_lanigera", "Chlorocebus_sabaeus", "Ciona_intestinalis", "Colobus_angolensis", "Cricetulus_griseus", "Cricetulus_griseus", "Danio_rerio", "Dasypus_novemcinctus", "Dipodomys_ordii", "Drosophila_melanogaster", "Eptatretus_burgeri", "Equus_caballus", "Felis_catus", "Ficedula_albicollis", "Fukomys_damarensis", "Gadus_morhua", "Gallus_gallus", "Gasterosteus_aculeatus", "Gorilla_gorilla", "Heterocephalus_glaber", "Homo_sapiens", "Ictidomys_tridecemlineatus", "Jaculus_jaculus", "Latimeria_chalumnae", "Lepisosteus_oculatus", "Loxodonta_africana", "Macaca_fascicularis", "Macaca_mulatta", "Macaca_nemestrina", "Mandrillus_leucophaeus", "Meleagris_gallopavo", "Mesocricetus_auratus", "Microcebus_murinus", "Microtus_ochrogaster", "Monodelphis_domestica", "Mus_caroli", "Mus_musculus", "Mus_pahari", "Mus_spretus", "Mustela_putorius", "Myotis_lucifugus", "Nannospalax_galili", "Nomascus_leucogenys", "Octodon_degus", "Ornithorhynchus_anatinus", "Oryctolagus_cuniculus", "Oryzias_latipes", "Otolemur_garnettii", "Ovis_aries", "Pan_paniscus", "Pan_troglodytes", "Panthera_pardus", "Panthera_tigris", "Papio_anubis", "Pelodiscus_sinensis", "Peromyscus_maniculatus", "Petromyzon_marinus", "Propithecus_coquereli", "Rattus_norvegicus", "Rhinopithecus_bieti", "Rhinopithecus_roxellana", "Saimiri_boliviensis", "Sarcophilus_harrisii", "Sus_scrofa", "Taeniopygia_guttata", "Takifugu_rubripes", "Tetraodon_nigroviridis", "Xenopus_tropicalis", "Xiphophorus_maculatus", "Acyrthosiphon_pisum", "Aedes_aegypti", "Amphimedon_queenslandica", "Anopheles_darlingi", "Anopheles_gambiae", "Apis_mellifera", "Atta_cephalotes", "Belgica_antarctica", "Bombyx_mori", "Brugia_malayi", "Caenorhabditis_brenneri", "Caenorhabditis_briggsae", "Caenorhabditis_japonica", "Caenorhabditis_remanei", "Capitella_teleta", "Crassostrea_gigas", "Culex_quinquefasciatus", "Danaus_plexippus", "Daphnia_pulex", "Dendroctonus_ponderosae", "Drosophila_ananassae", "Drosophila_erecta", "Drosophila_grimshawi", "Drosophila_mojavensis", "Drosophila_persimilis", "Drosophila_pseudoobscura", "Drosophila_sechellia", "Drosophila_simulans", "Drosophila_virilis", "Drosophila_willistoni", "Drosophila_yakuba", "Heliconius_melpomene", "Helobdella_robusta", "Lepeophtheirus_salmonis", "Lucilia_cuprina", "Mayetiola_destructor", "Melitaea_cinxia", "Mnemiopsis_leidyi", "Nematostella_vectensis", "Octopus_bimaculoides", "Pediculus_humanus", "Pristionchus_pacificus", "Rhodnius_prolixus", "Solenopsis_invicta", "Stegodyphus_mimosarum", "Strigamia_maritima", "Tribolium_castaneum", "Trichoplax_adhaerens", "Zootermopsis_nevadensis", "Acropora_digitifera", "Aedes_albopictus", "Apis_cerana", "Apis_dorsata", "Apis_florea", "Aplysia_californica", "Bactrocera_dorsalis", "Bactrocera_oleae", "Bemisia_tabaci", "Bicyclus_anynana", "Branchiostoma_belcheri", "Branchiostoma_floridae", "Camponotus_floridanus", "Ceratina_calcarata", "Ceratitis_capitata", "Ceratosolen_solmsi", "Cimex_lectularius", "Crassostrea_virginica", "Cyphomyrmex_costatus", "Drosophila_arizonae", "Drosophila_biarmipes", "Drosophila_bipectinata", "Drosophila_busckii", "Drosophila_elegans", "Drosophila_eugracilis", "Drosophila_ficusphila", "Drosophila_kikkawai", "Drosophila_miranda", "Drosophila_navojoa", "Drosophila_obscura", "Drosophila_rhopaloa", "Drosophila_serrata", "Drosophila_suzukii", "Drosophila_takahashii", "Dufourea_novaeangliae", "Echinococcus_granulosus", "Eufriesea_mexicana", "Eurytemora_affinis", "Metaseiulus_occidentalis", "Habropoda_laboriosa", "Hydra_vulgaris", "Limulus_polyphemus", "Linepithema_humile", "Mizuhopecten_yessoensis", "Monomorium_pharaonis", "Musca_domestica", "Myzus_persicae", "Nicrophorus_vespilloides", "Nilaparvata_lugens", "Onthophagus_taurus", "Montastraea_faveolata", "Orussus_abietinus", "Papilio_machaon", "Papilio_polytes", "Papilio_xuthus", "Parasteatoda_tepidariorum", "Pieris_rapae", "Pogonomyrmex_barbatus", "Priapulus_caudatus", "Pseudomyrmex_gracilis", "Saccoglossus_kowalevskii", "Stomoxys_calcitrans", "Trachymyrmex_cornetzi", "Trachymyrmex_septentrionalis", "Trachymyrmex_zeteki", "Varroa_destructor", "Vollenhovia_emeryi", "Wasmannia_auropunctata", "Bactrocera_cucurbitae", "Acinonyx_jubatus", "Balaenoptera_acutorostrata", "Bison_bison", "Bos_indicus", "Bos_mutus", "Bubalus_bubalis", "Camelus_bactrianus", "Camelus_dromedarius", "Camelus_ferus", "Castor_canadensis", "Ceratotherium_simum", "Chrysochloris_asiatica", "Condylura_cristata", "Echinops_telfairi", "Elephantulus_edwardii", "Eptesicus_fuscus", "Equus_asinus", "Equus_przewalskii", "Erinaceus_europaeus", "Hipposideros_armiger", "Leptonychotes_weddellii", "Lipotes_vexillifer", "Manis_javanica", "Marmota_marmota", "Meriones_unguiculatus", "Miniopterus_natalensis", "Myotis_brandtii", "Myotis_davidii", "Monachus_monachus", "Ochotona_princeps", "Odobenus_rosmarus", "Odocoileus_virginianus", "Orcinus_orca", "Orycteropus_afer", "Pantholops_hodgsonii", "Phascolarctos_cinereus", "Physeter_catodon", "Pongo_abelii", "Pteropus_alecto", "Pteropus_vampyrus", "Rhinolophus_sinicus", "Rousettus_aegyptiacus", "Sorex_araneus", "Trichechus_manatus", "Tursiops_truncatus", "Ursus_maritimus", "Vicugna_pacos", "Acanthisitta_chloris", "Acanthochromis_polyacanthus", "Alligator_mississippiensis", "Alligator_sinensis", "Anser_cygnoides", "Caprimulgus_carolinensis", "Apaloderma_vittatum", "Aptenodytes_forsteri", "Aquila_chrysaetos", "Balearica_regulorum", "Boleophthalmus_pectinirostris", "Buceros_rhinoceros", "Calidris_pugnax", "Callorhinchus_milii", "Calypte_anna", "Cariama_cristata", "Chaetura_pelagica", "Charadrius_vociferus", "Chelonia_mydas", "Chlamydotis_macqueenii", "Chrysemys_picta", "Clupea_harengus", "Colius_striatus", "Columba_livia", "Corvus_brachyrhynchos", "Coturnix_japonica", "Crocodylus_porosus", "Cuculus_canorus", "Cyprinodon_variegatus", "Cyprinus_carpio", "Egretta_garzetta", "Esox_lucius", "Eurypyga_helias", "Falco_cherrug", "Falco_peregrinus", "Fulmarus_glacialis", "Fundulus_heteroclitus", "Gavia_stellata", "Gavialis_gangeticus", "Gekko_japonicus", "Geospiza_fortis", "Haliaeetus_albicilla", "Haliaeetus_leucocephalus", "Haplochromis_burtoni", "Hippocampus_comes", "Ictalurus_punctatus", "Labrus_bergylta", "Larimichthys_crocea", "Lates_calcarifer", "Lepidothrix_coronata", "Leptosomus_discolor", "Lonchura_striata", "Manacus_vitellinus", "Maylandia_zebra", "Melopsittacus_undulatus", "Merops_nubicus", "Mesitornis_unicolor", "Nanorana_parkeri", "Neolamprologus_brichardi", "Nestor_notabilis", "Nipponia_nippon", "Nothobranchius_furzeri", "Notothenia_coriiceps", "Numida_meleagris", "Oncorhynchus_mykiss", "Opisthocomus_hoazin", "Parus_major", "Pelecanus_crispus", "Phaethon_lepturus", "Phalacrocorax_carbo", "Picoides_pubescens", "Poecilia_mexicana", "Pogona_vitticeps", "Protobothrops_mucrosquamatus", "Pseudopodoces_humilis", "Pterocles_gutturalis", "Pundamilia_nyererei", "Pygocentrus_nattereri", "Pygoscelis_adeliae", "Salmo_salar", "Scleropages_formosus", "Serinus_canaria", "Seriola_dumerili", "Sinocyclocheilus_anshuiensis", "Sinocyclocheilus_grahami", "Sinocyclocheilus_rhinocerous", "Stegastes_partitus", "Struthio_camelus", "Sturnus_vulgaris", "Tauraco_erythrolophus", "Thamnophis_sirtalis", "Tyto_alba", "Xenopus_laevis", "Zonotrichia_albicollis", "Ursus_arctos"])

totalDistance = 0
animal_tree = tree.copy()
### calculate total time over tree- losses are loss/total time. levelorder (default): every node on a level before is visited going to a lower level
for node in animal_tree.iter_descendants():
	totalDistance += node.dist	
# 	print(node.dist, totalDistance, len([leaf for leaf in node]))


###we start with the distribution of pfam metric values found in the youngest category of extant pfams (pfams < 500) 
filepath = '/Users/jennyjames/Desktop/LossRates/LossRatesForManuscript/PfamMLLikeBayesPipeline_MLLossAnimal_iteration4.txt'
Pfamfilepath = '/Users/jennyjames/Desktop/LossRates/LossRatesForManuscript/PfamUIDsTable_EnsemblAndNCBI.txt'

with open(filepath, 'r') as datafile, open(Pfamfilepath, 'r') as pfamdata:
	pfamdata = pfamdata.readlines()
# 	print(pfamdata[0])
	pfamdata = pfamdata[1:]
# 	print(len(pfamdata))
	pfamdata = [x.split('\t') for x in pfamdata]
	young = [x for x in pfamdata if float(x[8]) < 100]
	
	### we will exclude from our simulations pfams with no ISD or Clustering score
	### Animal ISD = 11, Animal Clustering = 15
	young = [x for x in young if x[11] != 'None' and x[15] != 'None']
	rnd_indices = np.random.choice(len(young), size=5000, replace = True)
	rnd_indices = list(rnd_indices)
	young_random_set = [young[x] for x in rnd_indices]
	
	
""" 
Regression slope results to use in simulations

ISD
breakpoint = 0.18

> slope(ISDBreak)
$MeanISD_AnimalSpecific
          Est. St.Err. t value CI(95%).l CI(95%).u
slope1 -3.4206 0.84015 -4.0714   -5.0675   -1.7736
slope2  2.8021 0.19105 14.6670    2.4276    3.1766

> intercept(ISDBreak)
$MeanISD_AnimalSpecific
              Est.
intercept1 -6.7068
intercept2 -7.8152


Clustering
breakpoint = 0.853

> slope(ClusteringBreak)
$MeanClustering_AnimalSpecific
          Est. St.Err. t value CI(95%).l CI(95%).u
slope1 -2.9891 0.34195 -8.7413   -3.6594   -2.3187
slope2  1.4237 0.13872 10.2630    1.1517    1.6956

> intercept(ClusteringBreak)
$MeanClustering_AnimalSpecific
              Est.
intercept1 -4.6896
intercept2 -8.4536



lambda = 0.007522682 
bc.backtransform <- function(y,L){(L*y+1)^(1/L)}
	 
"""	

pfam_bin_1 = []
pfam_bin_2 = []
pfam_bin_3 = []
pfam_bin_4 = []
pfam_bin_5 = []
pfam_bin_6 = []
pfam_bin_7 = []
pfam_bin_8 = []
pfam_bin_9 = []

### deterministically set loss rates, depending on ISD regression slopes
### loss rate is transformed, and the relationship is non-linear:
for pfam in young_random_set:

	if float(pfam[11]) > 0.1782:
		loss_rate_ISD = float(pfam[11]) * 2.8021 + -7.8152
	else:		
		loss_rate_ISD = float(pfam[11]) * -3.4206 + -6.7068
					
	###this is the backtransform for our data- lambda = 0.007522682 , lambda2 = 9.731337e-06
	###the unit is now loss rate/MY		
	loss_rate_ISD = (loss_rate_ISD * 0.007522682 + 1)**(1/0.007522682) - 9.731337e-06

	###How much time until a pfam is lost? Use a random draw from an exponential distribution, parameterised by loss rate/MY
	loss_rate_ran_time = np.random.exponential(scale=1/loss_rate_ISD, size=1)	
	
	###
	if loss_rate_ran_time > 10.0: 
		pfam_bin_1.append([float(pfam[11]), loss_rate_ISD])

	if loss_rate_ran_time > 50.0:
		pfam_bin_2.append([float(pfam[11]), loss_rate_ISD])

	if loss_rate_ran_time > 100.0:
		pfam_bin_3.append([float(pfam[11]), loss_rate_ISD])

	if loss_rate_ran_time > 250.0: 
		pfam_bin_4.append([float(pfam[11]), loss_rate_ISD])  
		
	if loss_rate_ran_time > 500.0: 
		pfam_bin_5.append([float(pfam[11]), loss_rate_ISD])
		
	if loss_rate_ran_time > 750.0: 
		pfam_bin_6.append([float(pfam[11]), loss_rate_ISD])
	
	if loss_rate_ran_time > 1000.0: 
		pfam_bin_7.append([float(pfam[11]), loss_rate_ISD])

	if loss_rate_ran_time > 1500.0: 
		pfam_bin_8.append([float(pfam[11]), loss_rate_ISD])

	if loss_rate_ran_time > 2000.0: 
		pfam_bin_9.append([float(pfam[11]), loss_rate_ISD])

print(len(pfam_bin_1))
print(len(pfam_bin_2))
print(len(pfam_bin_3))
print(len(pfam_bin_4))
print(len(pfam_bin_5))
print(len(pfam_bin_6))
print(len(pfam_bin_7))
print(len(pfam_bin_8))
print(len(pfam_bin_9))


print("c(10,"+ str(np.mean([float(x[0]) for x in pfam_bin_1])) +"," + str(np.std([float(x[0]) for x in pfam_bin_1])/np.sqrt(len(pfam_bin_1))) +"," + str(np.std([float(x[0]) for x in pfam_bin_1]))+"," + str(np.var([float(x[0]) for x in pfam_bin_1])) +")," )
print("c(50,"+ str(np.mean([float(x[0]) for x in pfam_bin_2])) +"," + str(np.std([float(x[0]) for x in pfam_bin_2])/np.sqrt(len(pfam_bin_2))) +"," + str(np.std([float(x[0]) for x in pfam_bin_2]))+"," + str(np.var([float(x[0]) for x in pfam_bin_2])) +")," )
print("c(100,"+ str(np.mean([float(x[0]) for x in pfam_bin_3])) +"," + str(np.std([float(x[0]) for x in pfam_bin_3])/np.sqrt(len(pfam_bin_3))) +"," + str(np.std([float(x[0]) for x in pfam_bin_3]))+"," + str(np.var([float(x[0]) for x in pfam_bin_3])) +")," )
print("c(250,"+ str(np.mean([float(x[0]) for x in pfam_bin_4])) +"," + str(np.std([float(x[0]) for x in pfam_bin_4])/np.sqrt(len(pfam_bin_4))) +"," + str(np.std([float(x[0]) for x in pfam_bin_4]))+"," + str(np.var([float(x[0]) for x in pfam_bin_4])) +")," )
print("c(500,"+ str(np.mean([float(x[0]) for x in pfam_bin_5])) +"," + str(np.std([float(x[0]) for x in pfam_bin_5])/np.sqrt(len(pfam_bin_5))) +"," + str(np.std([float(x[0]) for x in pfam_bin_5]))+"," + str(np.var([float(x[0]) for x in pfam_bin_5])) +")," )
print("c(750,"+ str(np.mean([float(x[0]) for x in pfam_bin_6])) +"," + str(np.std([float(x[0]) for x in pfam_bin_6])/np.sqrt(len(pfam_bin_6))) +"," + str(np.std([float(x[0]) for x in pfam_bin_6]))+"," + str(np.var([float(x[0]) for x in pfam_bin_6])) +"),")
print("c(1000,"+ str(np.mean([float(x[0]) for x in pfam_bin_7])) +"," + str(np.std([float(x[0]) for x in pfam_bin_7])/np.sqrt(len(pfam_bin_7))) +"," + str(np.std([float(x[0]) for x in pfam_bin_7]))+"," + str(np.var([float(x[0]) for x in pfam_bin_7])) +")," )
print("c(1500,"+ str(np.mean([float(x[0]) for x in pfam_bin_8])) +"," + str(np.std([float(x[0]) for x in pfam_bin_8])/np.sqrt(len(pfam_bin_8))) +"," + str(np.std([float(x[0]) for x in pfam_bin_8]))+"," + str(np.var([float(x[0]) for x in pfam_bin_8])) +")," )
print("c(2000,"+ str(np.mean([float(x[0]) for x in pfam_bin_9])) +"," + str(np.std([float(x[0]) for x in pfam_bin_9])/np.sqrt(len(pfam_bin_9))) +"," + str(np.std([float(x[0]) for x in pfam_bin_9]))+"," + str(np.var([float(x[0]) for x in pfam_bin_9])) +")" )

	
		
pfam_bin_1 = []
pfam_bin_2 = []
pfam_bin_3 = []
pfam_bin_4 = []
pfam_bin_5 = []
pfam_bin_6 = []
pfam_bin_7 = []
pfam_bin_8 = []
pfam_bin_9 = []

### deterministically set loss rates, depending on Clustering regression slopes
for pfam in young_random_set:

	if float(pfam[15]) > 0.853:
		loss_rate_ISD = float(pfam[15]) * 1.4237 + -8.4536
	else:		
		loss_rate_ISD = float(pfam[15]) * -2.9891 + -4.6896
					
	###this is the backtransform for our data- lambda = 0.007522682 , lambda2 = 9.731337e-06
	###the unit is now loss rate/MY		
	loss_rate_ISD = (loss_rate_ISD * 0.007522682 + 1)**(1/0.007522682) - 9.731337e-06

	###How much time until a pfam is lost? Use a random draw from an exponential distribution, parameterised by loss rate/MY
	loss_rate_ran_time = np.random.exponential(scale=1/loss_rate_ISD, size=1)	
	
	###
	if loss_rate_ran_time > 10.0: 
		pfam_bin_1.append([float(pfam[15]), loss_rate_ISD])

	if loss_rate_ran_time > 50.0:
		pfam_bin_2.append([float(pfam[15]), loss_rate_ISD])

	if loss_rate_ran_time > 100.0:
		pfam_bin_3.append([float(pfam[15]), loss_rate_ISD])

	if loss_rate_ran_time > 250.0: 
		pfam_bin_4.append([float(pfam[15]), loss_rate_ISD])  
		
	if loss_rate_ran_time > 500.0: 
		pfam_bin_5.append([float(pfam[15]), loss_rate_ISD])
		
	if loss_rate_ran_time > 750.0: 
		pfam_bin_6.append([float(pfam[15]), loss_rate_ISD])
	
	if loss_rate_ran_time > 1000.0: 
		pfam_bin_7.append([float(pfam[15]), loss_rate_ISD])

	if loss_rate_ran_time > 1500.0: 
		pfam_bin_8.append([float(pfam[15]), loss_rate_ISD])

	if loss_rate_ran_time > 2000.0: 
		pfam_bin_9.append([float(pfam[15]), loss_rate_ISD])

print(len(pfam_bin_1))
print(len(pfam_bin_2))
print(len(pfam_bin_3))
print(len(pfam_bin_4))
print(len(pfam_bin_5))
print(len(pfam_bin_6))
print(len(pfam_bin_7))
print(len(pfam_bin_8))
print(len(pfam_bin_9))


###standard error is calculated by dividing the standard deviation by the sample size's square root
print("c(10,"+ str(np.mean([float(x[0]) for x in pfam_bin_1])) +"," + str(np.std([float(x[0]) for x in pfam_bin_1])/np.sqrt(len(pfam_bin_1))) +"," + str(np.std([float(x[0]) for x in pfam_bin_1])) +"," + str(np.var([float(x[0]) for x in pfam_bin_1])) +")," )
print("c(50,"+ str(np.mean([float(x[0]) for x in pfam_bin_2])) +"," + str(np.std([float(x[0]) for x in pfam_bin_2])/np.sqrt(len(pfam_bin_2))) +"," + str(np.std([float(x[0]) for x in pfam_bin_2]))+"," + str(np.var([float(x[0]) for x in pfam_bin_2])) +")," )
print("c(100,"+ str(np.mean([float(x[0]) for x in pfam_bin_3])) +"," + str(np.std([float(x[0]) for x in pfam_bin_3])/np.sqrt(len(pfam_bin_3))) +"," + str(np.std([float(x[0]) for x in pfam_bin_3]))+"," + str(np.var([float(x[0]) for x in pfam_bin_3])) +")," )
print("c(250,"+ str(np.mean([float(x[0]) for x in pfam_bin_4])) +"," + str(np.std([float(x[0]) for x in pfam_bin_4])/np.sqrt(len(pfam_bin_4)))+"," + str(np.std([float(x[0]) for x in pfam_bin_4])) +"," + str(np.var([float(x[0]) for x in pfam_bin_4])) +")," )
print("c(500,"+ str(np.mean([float(x[0]) for x in pfam_bin_5])) +"," + str(np.std([float(x[0]) for x in pfam_bin_5])/np.sqrt(len(pfam_bin_5))) +"," + str(np.std([float(x[0]) for x in pfam_bin_5]))+"," + str(np.var([float(x[0]) for x in pfam_bin_5])) +")," )
print("c(750,"+ str(np.mean([float(x[0]) for x in pfam_bin_6])) +"," + str(np.std([float(x[0]) for x in pfam_bin_6])/np.sqrt(len(pfam_bin_6)))+"," + str(np.std([float(x[0]) for x in pfam_bin_6])) +"," + str(np.var([float(x[0]) for x in pfam_bin_6])) +"),")
print("c(1000,"+ str(np.mean([float(x[0]) for x in pfam_bin_7])) +"," + str(np.std([float(x[0]) for x in pfam_bin_7])/np.sqrt(len(pfam_bin_7))) +"," + str(np.std([float(x[0]) for x in pfam_bin_7]))+"," + str(np.var([float(x[0]) for x in pfam_bin_7])) +")," )
print("c(1500,"+ str(np.mean([float(x[0]) for x in pfam_bin_8])) +"," + str(np.std([float(x[0]) for x in pfam_bin_8])/np.sqrt(len(pfam_bin_8))) +"," + str(np.std([float(x[0]) for x in pfam_bin_8]))+"," + str(np.var([float(x[0]) for x in pfam_bin_8])) +")," )
print("c(2000,"+ str(np.mean([float(x[0]) for x in pfam_bin_9])) +"," + str(np.std([float(x[0]) for x in pfam_bin_9])/np.sqrt(len(pfam_bin_9))) +"," + str(np.std([float(x[0]) for x in pfam_bin_9]))+"," + str(np.var([float(x[0]) for x in pfam_bin_9])) +")" )

	


