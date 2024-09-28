#!/usr/bin/python
# -*- coding: UTF-8 -*-
import pyfastx

def peptideMotifScore(amino_file,codon_exp_file,output_file,occurrence_threshold=100):
    ###########################################################################################
    # 1_load amino fasta file
    ###########################################################################################
    amino = pyfastx.Fasta(amino_file)
    
    ###########################################################################################
    # 2_get density for tripepetide density
    ###########################################################################################
    tripeptide_seq =  {}

    with open(codon_exp_file,'r') as input:
        for line in input:
            trans_id,codon_pos,exp = line.split()
           
            # check transid in AA fasta file (revised code)!
            if trans_id in amino.keys():
                # 3nt window to fetch tri-peptide motif
                
                seq = str(amino[trans_id])
            
                if int(codon_pos) - 1 >= 1 and int(codon_pos) + 1 <= len(seq):
                    type1 = amino.fetch(trans_id,(int(codon_pos) - 1,int(codon_pos) + 1))
                else:
                    type1 = ""
                
                if int(codon_pos) - 2 >= 1 and int(codon_pos) <= len(seq):
                    type2 = amino.fetch(trans_id,(int(codon_pos) - 2,int(codon_pos)))
                else:
                    type2 = ""
                    
                if int(codon_pos) >= 1 and int(codon_pos) + 2 <= len(seq):
                    type3 = amino.fetch(trans_id,(int(codon_pos),int(codon_pos) + 2))
                else:
                    type3 = ""
                
                tripeptide_motif_list = [type1, type2, type3]
                
                # filter motif length is 3
                tripeptide_motif_list_filterd = [i for i in tripeptide_motif_list if len(i) == 3]
                
                # get tripeptide and save in dict
                for motif in tripeptide_motif_list_filterd:
                    if motif in tripeptide_seq:
                        tmp_exp,tmp_count = tripeptide_seq[motif]
                        tripeptide_seq[motif] = [float(exp) + tmp_exp,tmp_count + 1]
                    else:
                        tripeptide_seq[motif] = [float(exp),1]
                        
    ###########################################################################################
    # 3_filter motif occurrence
    ###########################################################################################
    tripeptide_filtered = {}

    # filter motifs with more than 100 occurrences
    for key,val in tripeptide_seq.items():
        if val[1] > occurrence_threshold:
            tripeptide_filtered[key] = [val[0],val[1]]
            
    ###########################################################################################
    # 4_output
    ###########################################################################################
    total_motifs = 0
    total_exp = 0
    for key,val in tripeptide_filtered.items():
        total_motifs += val[1]
        total_exp += val[0]
        
    # output
    out_file = open(output_file, 'w',encoding='utf-8')
    for key,val in tripeptide_filtered.items():
        # norm_val = (val[0]/(val[1]/total_motifs))/1000000
        norm_val = (val[0]/total_exp)*1000
        # norm_val = (val[0]/total_exp)*(val[1]/total_motifs)
        out_file.write("\t".join([str(key),str(norm_val)]) + '\n')
    out_file.close()