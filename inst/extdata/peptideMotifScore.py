#!/usr/bin/python
# -*- coding: UTF-8 -*-
import pyfastx

def peptideMotifScore(amino_file,codon_exp_file,output_file,occurrence_threshold):
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
            
            # # get tripeptide position
            # if int(codon_pos)%3 == 0:
            #     pos = int(codon_pos)//3
            # else:
            #     pos = int(int(codon_pos)//3) + 1
                
            # # check transid in AA fasta file
            # if trans_id in amino.keys():
            #     tripeptide_motif = amino.fetch(trans_id,(pos*3 - 2,pos*3))
            
            #     # make sure length is 3
            #     if len(tripeptide_motif) == 3:
            #         # get tripeptide and save in dict
            #         if tripeptide_motif in tripeptide_seq:
            #             tmp_exp,tmp_count = tripeptide_seq[tripeptide_motif]
            #             tripeptide_seq[tripeptide_motif] = [float(exp) + tmp_exp,tmp_count + 1]
            #         else:
            #             tripeptide_seq[tripeptide_motif] = [float(exp),1]
                
            # check transid in AA fasta file (revised code)!
            if trans_id in amino.keys():
                # 3nt window to fetch tri-peptide motif
                
                # type1 = amino.fetch(trans_id,(int(codon_pos) - 1,int(codon_pos) + 1))
                # type2 = amino.fetch(trans_id,(int(codon_pos) - 2,int(codon_pos)))
                # type3 = amino.fetch(trans_id,(int(codon_pos),int(codon_pos) + 2))
                
                try:
                    type1 = amino.fetch(trans_id,(int(codon_pos) - 1,int(codon_pos) + 1))
                except Exception as e:
                    continue
                    
                try:
                    type2 = amino.fetch(trans_id,(int(codon_pos) - 2,int(codon_pos)))
                except Exception as e:
                    continue
                
                try:
                    type3 = amino.fetch(trans_id,(int(codon_pos),int(codon_pos) + 2))
                except Exception as e:
                    continue
                
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
    for key,val in tripeptide_filtered.items():
        total_motifs += val[1]
        
    # output
    out_file = open(output_file, 'w',encoding='utf-8')
    for key,val in tripeptide_filtered.items():
        norm_val = (val[0]/(val[1]/total_motifs))/1000000
        out_file.write("\t".join([str(key),str(norm_val)]) + '\n')
    out_file.close()