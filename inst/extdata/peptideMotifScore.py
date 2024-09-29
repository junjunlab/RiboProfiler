import pyfastx

def peptideMotifScore(amino_file,longest_trans_file,normed_file,min_counts,average_normalization,output_file,occurrence_threshold=100):
    ###########################################################################################
    # 1_load amino acids
    ###########################################################################################
    amino = pyfastx.Fasta(amino_file)
    
    ###########################################################################################
    # 2_construct whole gene pos
    ###########################################################################################
    gene_info = {}
    exp_whole_trans = {}

    with open(longest_trans_file,'r') as input:
        for line in input:
            _, _, _, transid, _, _, _, _, utr5, cds, utr3  = line.split()
            
            if int(cds)%3 == 0:
                gene_info[transid] = [int(utr5),int(cds),int(utr3)]
                
                pos = [0 for i in range(1,int(int(cds)/3) + 1)]
                exp_whole_trans[transid] = pos
                
    ###########################################################################################
    # 3_save density
    ###########################################################################################
    gene_counts = {}
    with open(normed_file,'r') as input:
        for line in input:
            _,_,_,_,_,_,trans_pos,trans_id,counts,_,_,_,exp = line.split()
            utr5,cds,_ = gene_info[trans_id]
            cds_pos = int(trans_pos) - utr5
            # transform into codon pos
            if cds_pos >= 1 and cds_pos <= cds:
                if cds_pos%3 == 0:
                    codon_pos = int(cds_pos/3)
                else:
                    codon_pos = cds_pos//3 + 1
                        
                exp_whole_trans[trans_id][(codon_pos - 1)] = float(exp)
                
                # save gene counts
                if trans_id in gene_counts:
                    gene_counts[trans_id] += int(counts)
                else:
                    gene_counts[trans_id] = int(counts)
                
    # filter gene counts
    gene_counts_filtered = {}

    for key,val in gene_counts.items():
        if val > min_counts:
            gene_counts_filtered[key] = val
            
    ###########################################################################################
    # 4_average normaliation
    ###########################################################################################
    average_means = {}
    
    if average_normalization == "yes":
        for key,val in exp_whole_trans.items():
            if key in gene_counts_filtered:
                total_density = sum(val)
            
                if total_density > 0:
                    average_exp = total_density/len(val)
                    normed_val = [i/average_exp for i in val]
                else:
                    normed_val = normed_val
                    
                # save in new dict
                average_means[key] = normed_val
    else:
        for key,val in exp_whole_trans.items():
            if key in gene_counts_filtered:
                average_means[key] = val
        
    ###########################################################################################
    # 5_conut motif density along transcript
    ###########################################################################################
    aa_id = [id for id in amino.keys()]

    motif_occurance = {}

    for tid in average_means.keys():
        if tid in aa_id:
            # get amino acid
            sequence = str(amino[tid])
            combinations = [sequence[i:i+3] for i in range(len(sequence) - 2)]
            
            # get density
            density = exp_whole_trans[tid]
            density_motif = [sum(density[i:i+3]) for i in range(len(density) - 2)]
            
            # loop save motif and density
            for i in range(len(combinations)):
                motif = combinations[i]
                exp = density_motif[i]
                
                # save in dict
                if motif in motif_occurance:
                    motif_exp,count = motif_occurance[motif]
                    motif_occurance[motif] = [motif_exp + exp,count + 1]
                else:
                    motif_occurance[motif] = [exp,1]
                    
    ###########################################################################################
    # 6_output
    ###########################################################################################
    tripeptide_filtered = {}

    # filter motifs with more than 100 occurrences
    for key,val in motif_occurance.items():
        if val[1] > occurrence_threshold:
            tripeptide_filtered[key] = [val[0],val[1]]

    out_file = open(output_file, 'w')
    for key,val in tripeptide_filtered.items():
        norm_val = val[0]/val[1]
        out_file.write("\t".join([str(key),str(norm_val)]) + '\n')
    out_file.close()
        