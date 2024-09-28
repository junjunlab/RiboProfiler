import pyfastx

def averageMotifPauseScore(amino_file,longest_trans_file,input_file,output_file,norm_type,window=50,occurrence_threshold=100):
    ###########################################################################################
    # 1_load amino acids
    ###########################################################################################
    amino = pyfastx.Fasta(amino_file)
    
    ###########################################################################################
    # 2_load longtest transcript anno file
    ###########################################################################################
    gene_info = {}
    exp_whole_trans = {}

    with open(longest_trans_file,'r') as input:
        for line in input:
            _, _, _, transid, _, _, _, _, utr5, cds, utr3  = line.split()
            gene_info[transid] = [int(utr5),int(cds),int(utr3)]
            
            pos = [0 for i in range(1,int(cds) + 1)]
            exp_whole_trans[transid] = pos
    
    ###########################################################################################
    # 3_save for each transcript across every position
    ###########################################################################################
    with open(input_file,'r') as input:
        for line in input:
            _,_,_,_,_,_,trans_pos,trans_id,counts,_,_,_,exp = line.split()
            utr5,cds,_ = gene_info[trans_id]
            cds_pos = int(trans_pos) - utr5
            
            if cds_pos >= 1 and cds_pos <= cds:
                if norm_type == "rpm":
                    exp_whole_trans[trans_id][(cds_pos - 1)] = float(counts)
                else:
                    exp_whole_trans[trans_id][(cds_pos - 1)] = float(exp)

    ###########################################################################################
    # 4_get motif pause score in window
    ###########################################################################################
    pause_score = {}
    with open(input_file,'r') as input:
        for line in input:
            _,_,_,_,_,_,trans_pos,trans_id,counts,_,_,_,exp = line.split()
            utr5,cds,_ = gene_info[trans_id]
            cds_pos = int(trans_pos) - utr5
            
            # exclude first and last 50nt
            if cds_pos >= window + 1 and cds_pos <= cds - window:
                nt_exp = exp_whole_trans[trans_id]
                window_exp = nt_exp[(cds_pos - window - 1):(cds_pos + window)]
                
                # normalize average reads in window
                average_mean = sum(window_exp)/(window + window + 1)
                normalized_window_exp = [i/average_mean for i in window_exp]
                
                # get EPA tri amino acid score
                average_sum_epa = sum(normalized_window_exp[(window - 3 - 1):(window + 5)])
                
                # calculate codon position and get tri-amino acid motif
                if cds_pos%3 == 0:
                    codon_pos = int(cds_pos/3)
                else:
                    codon_pos = cds_pos//3 + 1
                    
                motif = amino.fetch(trans_id,(codon_pos - 1,codon_pos + 1))
                
                # save in dict
                if motif in pause_score:
                    tmp_exp,tmp_count = pause_score[motif]
                    pause_score[motif] = [average_sum_epa + tmp_exp,tmp_count + 1]
                else:
                    pause_score[motif] = [average_sum_epa,1]
                        
    ###########################################################################################
    # 5_filter motif occurrence
    ###########################################################################################
    tripeptide_filtered = {}

    # filter motifs with more than 100 occurrences
    for key,val in pause_score.items():
        if val[1] > occurrence_threshold:
            tripeptide_filtered[key] = [val[0],val[1]]
            
    ###########################################################################################
    # 6_output
    ###########################################################################################
    out_file = open(output_file, 'w')
    for key,val in tripeptide_filtered.items():
        # norm_val = val[0]
        norm_val = val[0]/val[1]
        out_file.write("\t".join([str(key),str(norm_val)]) + '\n')
    out_file.close()