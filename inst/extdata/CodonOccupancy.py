def codonOccupancy(cds_fasta_file,codon_pos_exp_file,output_file):
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
    # 1_prepare 64 codons Dict
    ######################################################################################################################
    import pyfastx

    # load cds fasta
    fa = pyfastx.Fasta(cds_fasta_file)

    codon_exp_dict = {}

    # load codon density file
    with open(codon_pos_exp_file,'r') as input:
        for line in input:
                trans_id,codon_pos,exp = line.split()
                codon_seq = fa.fetch(trans_id,(int(codon_pos)*3 - 2,int(codon_pos)*3))
                
                if codon_seq in codon_exp_dict:
                    tmp_exp,tmp_count = codon_exp_dict[codon_seq]
                    codon_exp_dict[codon_seq] = [float(exp) + tmp_exp,tmp_count + 1]
                else:
                    codon_exp_dict[codon_seq] = [float(exp),1]
                
    # get total codons
    total_codons = 0
    total_exp = 0
    for key,val in codon_exp_dict.items():
        total_codons += val[1]
        total_exp += val[0]
    
    # output
    out_file = open(output_file, 'w')
    for key,val in codon_exp_dict.items():
        # norm_val = (val[0]/(val[1]/total_codons))/1000000
        norm_val = (val[0]/total_exp)*1000
        anno = codonDict[key]
        out_file.write("\t".join([key,anno[0],anno[1],str(norm_val)]) + '\n')