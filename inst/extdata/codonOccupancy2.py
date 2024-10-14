import pyfastx
from collections import Counter

def codonOccupancy2(cds_fasta_file,codon_pos_exp_file,output_file):
    ######################################################################################################################
    # 1_prepare 64 codons Dict
    ######################################################################################################################
    baseATGC = ["A","T","G","C"]

    codonDict = {}

    for i in baseATGC:
        for j in baseATGC:
            for z in baseATGC:
                key = i+j+z
                if key in ["GGT", "GGC", "GGA", "GGG"]:
                    codonDict[key] = ["Gly","G"]
                elif key in ["GCC", "GCA", "GCG", "GCT"]:
                    codonDict[key] = ["Ala","A"]
                elif key in ["GTT", "GTC", "GTA", "GTG"]:
                    codonDict[key] = ["Val","V"]
                elif key in ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"]:
                    codonDict[key] = ["Leu","L"]
                elif key in ["ATT","ATC","ATA"]:
                    codonDict[key] = ["Ile","I"]
                elif key in ["CCT", "CCC", "CCA", "CCG"]:
                    codonDict[key] = ["Pro","P"]
                elif key in ["TTT", "TTC"]:
                    codonDict[key] = ["Phe","F"]
                elif key in ["TAT", "TAC"]:
                    codonDict[key] = ["Tyr","Y"]
                elif key in ["TGG"]:
                    codonDict[key] = ["Trp","W"]
                elif key in ["TCT", "TCA", "TCC", "TCG", "AGT", "AGC"]:
                    codonDict[key] = ["Ser","S"]
                elif key in ["ACT", "ACC", "ACG", "ACA"]:
                    codonDict[key] = ["Thr","T"]
                elif key in ["TGT", "TGC"]:
                    codonDict[key] = ["Cys","C"]
                elif key in ["ATG"]:
                    codonDict[key] = ["Met","M"]
                elif key in ["AAT","AAC"]:
                    codonDict[key] = ["Asn","N"]
                elif key in ["CAA","CAG"]:
                    codonDict[key] = ["Gln","Q"]
                elif key in ["GAT","GAC"]:
                    codonDict[key] = ["Asp","D"]
                elif key in ["GAA","GAG"]:
                    codonDict[key] = ["Glu","E"]
                elif key in ["AAA","AAG"]:
                    codonDict[key] = ["Lys","K"]
                elif key in ["CGT","CGC","CGG","CGA","AGA","AGG"]:
                    codonDict[key] = ["Arg","R"]
                elif key in ["CAT","CAC"]:
                    codonDict[key] = ["His","H"]
                elif key in ["TGA","TAG","TAA"]:
                    codonDict[key] = ["Stop","SP"]
               
    ######################################################################################################################
    # 2_calculate codon frequence
    ######################################################################################################################     
    # load cds fasta
    fa = pyfastx.Fastx(cds_fasta_file)
    
    codon_freq_dict = {}

    for tid,seq in fa:
        seq_len = len(seq)
        # filter cds length divided by 3
        if len(seq) % 3 == 0:
            codon_seq_list = []
            for i in range(1,seq_len,3):
                codon_seq_list.append(seq[(i-1):(i+2)])
                
            # count freq for codon
            codon_numbers = Counter(codon_seq_list)
            total_codons = len(codon_seq_list)
            
            # save codon frequence
            for codon,nm in codon_numbers.items():
                key = "|".join([tid,codon])
                codon_freq = nm/total_codons
                codon_freq_dict[key] = codon_freq
                
    ######################################################################################################################
    # 3_genen total expression
    ######################################################################################################################  
    gene_total_exp = {}

    # load codon density file
    with open(codon_pos_exp_file,'r') as input:
        for line in input:
            trans_id,_,exp = line.split()
            
            if trans_id in gene_total_exp:
                gene_total_exp[trans_id] += float(exp)
            else:
                gene_total_exp[trans_id] = float(exp)
            
    ######################################################################################################################
    # 4_calculate codon dewell time
    ###################################################################################################################### 
    # load cds fasta
    cdsfa = pyfastx.Fasta(cds_fasta_file)

    codon_rod_dict = {}

    # load codon density file
    with open(codon_pos_exp_file,'r') as input:
        for line in input:
            trans_id,codon_pos,exp = line.split()
            codon_seq = cdsfa.fetch(trans_id,(int(codon_pos)*3 - 2,int(codon_pos)*3))
            
            key = "|".join([trans_id,codon_seq])
            codon_freq = codon_freq_dict[key]
            norm_exp = float(exp)/gene_total_exp[trans_id]
            
            rod = norm_exp/codon_freq
            
            # save
            if codon_seq in codon_rod_dict:
                tmp_rod,count = codon_rod_dict[codon_seq]
                codon_rod_dict[codon_seq] = [tmp_rod + rod,count + 1]
            else:
                codon_rod_dict[codon_seq] = [rod,1]
                
    ######################################################################################################################
    # 5_output
    ###################################################################################################################### 
    # output
    out_file = open(output_file, 'w')
    for key,val in codon_rod_dict.items():
        mean_val = val[0]/val[1]
        anno = codonDict[key]
        out_file.write("\t".join([key,anno[0],anno[1],str(mean_val)]) + '\n')
    out_file.close()