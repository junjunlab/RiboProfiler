import re
import pyfastx

def relDistTripeptideMotif(amino_file,longest_trans_file,normed_file,motif,output_file,upstream = -50,downstream = 50):
    ###########################################################################################
    # 1_read amino acids file
    ###########################################################################################
    amino = pyfastx.Fasta(amino_file)
    
    ###########################################################################################
    # 2_load longtest transcript anno file
    ###########################################################################################
    gene_info = {}

    with open(longest_trans_file,'r') as input:
        for line in input:
            _, _, _, transid, _, _, _, _, utr5, cds, utr3  = line.split()
            gene_info[transid] = [int(utr5),int(cds),int(utr3)]
            
    # create relative distance vector
    rel_dist= {}
    for i in range(int(upstream),int(downstream + 1)):
        rel_dist[i] = [0,0]
        
    ###########################################################################################
    # 2_load longtest transcript anno file
    ###########################################################################################
    with open(normed_file,'r') as input:
        for line in input:
            _,_,_,_,_,_,trans_pos,trans_id,counts,_,_,exp = line.split()
            
            if trans_id in gene_info:
                seq = str(amino[trans_id])
                
                # search motif patterns for amino acids
                pattern = re.compile(motif, flags=re.IGNORECASE)
                pos_res = pattern.finditer(seq)
                
                # motif codon position
                codon_motif_pos = [pos.start() + 1 for pos in pos_res]
                
                # check whether has target motif
                if len(codon_motif_pos) > 0:
                    # motif codon position to nt position
                    motif_nt_pos = [pos*3- 2 for pos in codon_motif_pos]
                    
                    # calculate relative distance to motif
                    utr5_len = gene_info[trans_id][0]
                    rel_pos_to_motif = [int(trans_pos) - utr5_len - nt_pos for nt_pos in motif_nt_pos]
                    
                    # filter rel_pos in limited range (-50,50)
                    rel_pos_to_motif_in_range = [i for i in rel_pos_to_motif if i >= int(upstream) and i <= int(downstream)]
                    
                    # save in dict
                    for j in rel_pos_to_motif_in_range:
                        if j in rel_dist:
                            rel_dist[j][0] += 1
                            rel_dist[j][1] += float(exp)
                            # rel_dist[j][1] += int(counts)
                        
    # normalize density | Average ribosome occupancy
    sum_counts = sum([val[1] for key,val in rel_dist.items()])
    average_exp = sum_counts/len(range(int(upstream),int(downstream + 1)))
    rel_dist_norm = {key:val[1]/average_exp for key,val in rel_dist.items()}
    
    ###########################################################################################
    # 3_output file
    ###########################################################################################
    # output
    out_file = open(output_file, 'w')

    for key,val in rel_dist_norm.items():
        out_file.write("\t".join([str(key),str(val)]) + '\n')