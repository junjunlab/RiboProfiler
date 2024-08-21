from pyfaidx import Fasta
from Bio.Seq import Seq

def fetchSequneceFromGeneAnno(gene_file,genome_file,output_file,type = "cds",coding_type = "NT",table = 1):
    genome = Fasta(genome_file)

    seq_Dict = {}
    aa_Dict = {}
   
    ######################################################################################################                
    # extract sequence
    ######################################################################################################
    with open(gene_file) as ganno:
            for line in ganno:
                # split lines
                _,_,_,tid,chrom,strand,cds_rg,exon_rg,utr5_len,cds_len,utr3_len = line.split("\t")
                
                # check seq type
                if type == "cds":
                    key = "".join([">",tid," ",cds_len])
                    rg_slice = cds_rg.split(",")
                else:
                    rg_slice = exon_rg.split(",")
                    key = "".join([">",tid," ",str(int(cds_len) + int(utr5_len) + int(utr3_len))])
                    
                seq_Dict[key] = ""
                
                lensplice = len(rg_slice)
                
                # incude stop codon
                if strand == "+":
                    st = rg_slice[lensplice - 1].split(":")[0]
                    sp = str(int(rg_slice[lensplice - 1].split(":")[1]) + 3)
                    rg_slice[lensplice - 1] = ":".join([st,sp])
                else:
                    st = str(int(rg_slice[0].split(":")[0]) - 3)
                    sp = rg_slice[lensplice - 1].split(":")[1]
                    rg_slice[0] = ":".join([st,sp])
                
                # extract seqeunce from genome
                for i in rg_slice:
                    seq = genome[chrom][(int(i.split(":")[0]) - 1):int(i.split(":")[1])]
                    seq_Dict[key] += str(seq)
                    
                # reverse_complement for negtive strand gene
                if strand == "-":
                    seq_negstrand = Seq(seq_Dict[key])
                    seq_Dict[key] = str(seq_negstrand.reverse_complement())
                    
                # translate into anmino acids
                if type == "cds" and coding_type == "AA":
                    # check length whether can be divided with 3
                    if len(seq_Dict[key]) % 3 == 0:
                        seq_coding = Seq(seq_Dict[key])
                        aa_seq = str(seq_coding.translate(table = int(table),to_stop = True))
                        aa_Dict[" ".join([">",tid," ",str(len(aa_seq))])] = aa_seq
    
    ######################################################################################################                
    # output file
    ######################################################################################################
    output_fa = open(output_file,'w')

    # separate sequences
    if type == "cds" and coding_type == "AA":
        final_output_dict = aa_Dict
    else:
        final_output_dict = seq_Dict
        
    # get seq
    for key,val in final_output_dict.items():
        if len(val) != 0:
            output_fa.write(key + '\n')
            while len(val) > 80:
                output_fa.write(val[0:80] + '\n')
                val = val[80:len(val)]
            output_fa.write(val + '\n')

    # file close
    output_fa.close()