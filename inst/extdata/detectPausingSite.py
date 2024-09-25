def detectPausingSite(longest_trans_file,normed_file,output_file,min_counts=64,window=500):
    ###########################################################################################
    # 1_load annotatin file
    ###########################################################################################
    gene_info = {}

    with open(longest_trans_file,'r') as input:
        for line in input:
            _, gene_name, _, transid, _, _, _, _, utr5, cds, utr3  = line.split()
            gene_info[transid] = [int(utr5),int(cds),gene_name]
            
    ###########################################################################################
    # 2_filter gene with low counts
    ###########################################################################################
    total_exp_filter = {}
    
    with open(normed_file,'r') as input:
        for line in input:
            # fileds = line.split()
            # print(line.split())
            _,_,_,_,_,_,trans_pos,trans_id,counts,_,_,exp = line.split()
            
            utr5,cds,_ = gene_info[trans_id]
            
            if int(trans_pos) >= int(utr5) and int(trans_pos) <= int(utr5) + int(cds):
                if trans_id in total_exp_filter:
                    total_exp_filter[trans_id] += int(counts)
                else:
                    total_exp_filter[trans_id] = 0
            
            
    # filter low conuts gene
    total_exp_retained = {key:val for key,val in total_exp_filter.items() if val >= int(min_counts)}
    
    ###########################################################################################
    # 3_load normed file and save in dict
    ###########################################################################################
    trans_id_dict = {}

    with open(normed_file,'r') as input:
        for line in input:
            _,_,_,_,_,_,trans_pos,trans_id,counts,_,_,exp = line.split()
            if trans_id in total_exp_retained:
                if trans_id in trans_id_dict:
                    # trans_id_dict[trans_id].update({int(trans_pos): float(exp)})
                    if int(trans_pos) in trans_id_dict[trans_id]:
                        trans_id_dict[trans_id][int(trans_pos)] += float(exp)
                    else:
                        trans_id_dict[trans_id][int(trans_pos)] = float(exp)
                else:
                    trans_id_dict[trans_id] = {int(trans_pos): float(exp)}
                    
    # sorted by transpos
    trans_id_sorted = {}
    for key,pos_dict in trans_id_dict.items():
        # sort
        sorted_keys = sorted(pos_dict)
        sorted_dict = dict(zip(sorted_keys, [pos_dict[key] for key in sorted_keys]))
        trans_id_sorted[key] = sorted_dict
        
    ###########################################################################################
    # 4_calculate pausing score according to pausepred methods
    ###########################################################################################
    for key,pos_dict in trans_id_sorted.items():
        pos_dict_new = {}
        
        # loop to calculate pausing score for each position
        for pos,exp in pos_dict.items():
            pos_n = pos + int(window)
            pos_n_mid = pos + int(window)/2
            pos_n_extend = pos + 1.5*int(window)
            
            sum_1 = sum([pp for pi,pp in pos_dict.items() if pi >= pos and pi <= pos_n])
            sum_2 = sum([pp for pi,pp in pos_dict.items() if pi >= pos_n_mid and pi <= pos_n_extend])
            
            if sum_1*sum_2 == 0:
                si = 0
            else:
                si = (int(window)/2)*exp*((sum_1 + sum_2)/(sum_1*sum_2))
            
            # save in dict
            pos_dict_new[pos] = [exp,si]
        
        # save in dict
        trans_id_sorted[key] = pos_dict_new
        
    ###########################################################################################
    # 5_output data
    ###########################################################################################
    # output
    out_file = open(output_file, 'w')

    for key,val in trans_id_sorted.items():
        for pos,pause in val.items():
            out_file.write("\t".join([str(gene_info[key][2]),str(key),str(pos),str(pause[0]),str(pause[1])]) + '\n')
    out_file.close()
            
