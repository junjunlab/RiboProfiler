import pyfastx

def codonRelativeDist(cds_fasta_file,codon_pos_exp_file,output_file):
    ######################################################################################################################
    # 1_save codon seq and relative dist
    ######################################################################################################################
    # load cds fasta
    fa = pyfastx.Fastx(cds_fasta_file)
    
    codon_dist = {}

    for tid,seq in fa:
        seq_len = len(seq)
        # filter cds length divided by 3
        if len(seq) % 3 == 0:
            codons = [seq[(i-1):(i+2)] for i in range(1, seq_len, 3)]
            len_codons = len(codons)
            rel_pos = [(codons[i-1],round(i/len_codons,2),0) for i in range(1,len_codons + 1)]
            # save
            codon_dist[tid] = rel_pos
            
    ######################################################################################################################
    # 2_calculate codon frequence
    ######################################################################################################################
   
    # load codon density file
    with open(codon_pos_exp_file,'r') as input:
        for line in input:
            trans_id,codon_pos,exp = line.split()
            
            # get codon information
            codon_seq,rel_pos,occup = codon_dist[trans_id][int(codon_pos) - 1]
            codon_dist[trans_id][int(codon_pos) - 1] = (codon_seq,rel_pos,float(exp))
            
    ######################################################################################################################
    # 3_output
    ###################################################################################################################### 
    out_file = open(output_file, 'w')

    for key,val in codon_dist.items():
        for i in val:
            out_file.write("\t".join([key,i[0],str(i[1]),str(i[2])]) + '\n')
    out_file.close()