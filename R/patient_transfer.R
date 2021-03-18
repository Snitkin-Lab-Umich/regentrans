#patient flow network 
#network and *igraph packages 
#cols: source_facil, dest_facil, n_transfers
#one of each pair, ask them to make sure that is true 
#often part of data is masked (anything under 10 = -9, assume it means 1) 
#user should change that, int > 0 (unless they want to weight, warning)

#to do: 
#function takes in edge list and snv distances facility of each isolate, optional patient
#user choose snv threshold (recommend multiple, see how it changes)
#outputs for each closely related pair of two isolates from different faiclities and different patients 
#return pair, snv distance, facility of each isolate, number of patient transfers between those two facilities 
#for faiclity pairs with closely related isolates, do they have a lot of patient transfers between them? 
#relationship between n closely related pairs and n patient transfers 
#directed? 

#in the get_snv_distances give optional thing, add the two cols on 
#showcase how to use them in the vignette 

#option of direct or indirect transfers 

#let people access those functions too 

#another function to summarize 
#summarize number of closely related paris for each facility 
#facil_1, Facil_2, transfer_info, n_closely_related_pairs

#transfer is just a number 

#goal: for certain set of closely related isolates from different facilities, 
#we want to get the number of direct transfers between those facilities 

#edge matrix will have more facilites than you have isolates from 
#include that in a test
#only include facilites we have samples from 
#can plot scatterploy of transfer info vs. n closely related paris 

