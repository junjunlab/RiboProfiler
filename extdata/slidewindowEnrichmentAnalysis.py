# slidewindow function to smooth density ratio
def SlideWindowEnrich(inputDict):
    # calculate window +-20
    sum1 = 0
    sum2 = 0
    start = 1
    end = len(inputDict)
    # start_sum = 21
    start_sum = start + 20
    # end_sum = 4578139
    end_sum = end - 20

    outputDict = {}
    
    # slide window calculate
    # lines 1-21 sum 1-(J+20)
    for J in range(start, start_sum + 1):
        for X in range(start, J+20+1):
            sum1 += float(inputDict[X][0])
            sum2 += float(inputDict[X][1])
        if sum1 != 0:
            ratio = sum2 / sum1  
        else:
            ratio = 0.0
        outputDict[J] = ratio
        sum1 = 0
        sum2 = 0

    # lines 22-(end-19) sum (K-20)-(K+20)
    for K in range(start_sum + 1, end_sum + 1):
        for Y in range(K-20, K+20+1):
            sum1 += float(inputDict[Y][0])
            sum2 += float(inputDict[Y][1])
        if sum1 != 0:
            ratio = sum2 / sum1   
        else:
            ratio = 0.0
        outputDict[K] = ratio
        sum1 = 0
        sum2 = 0

    # lines 4578140-end
    for L in range(end_sum + 1, end + 1):
        for Z in range(L-20, end + 1):
            sum1 += float(inputDict[Z][0])
            sum2 += float(inputDict[Z][1])
        if sum1 != 0:
            ratio = sum2 / sum1
        else:
            ratio = 0.0
        outputDict[L] = ratio
        sum1 = 0
        sum2 = 0
    return(outputDict)

# calculate enrichmentvalue
def SingleGeneEnrichment(inputIPFile,inputInputFile,geneName):
    inter_Dict = {}
    trans_Dict = {}
    inter_gene_Dict = {}
    trans_gene_Dict = {}
    combined_Dict = {}

    # 1.open denisty file and save in dict
    with open(inputIPFile,'r') as interactome, \
        open(inputInputFile,'r') as translatome:
        # opne IP file
        for line in interactome:
            fileds = line.split()
            gene_name = fileds[0].split('|')[0]
            if gene_name == geneName:
                key = ':'.join([fileds[0],str(fileds[1])])
                inter_Dict[key] = float(fileds[3])
                # prepare gene cds info
                length = int(fileds[0].split('|')[4]) + int(fileds[0].split('|')[5]) + 2
                inter_gene_Dict = {':'.join([fileds[0],str(i)]):'' for i in range(1,length + 1)}

        # opne Input file
        for line in translatome:
            fileds = line.split()
            gene_name = fileds[0].split('|')[0]
            if gene_name == geneName:
                key = ':'.join([fileds[0],str(fileds[1])])
                trans_Dict[key] = float(fileds[3])
                # prepare gene cds info
                length = int(fileds[0].split('|')[4]) + int(fileds[0].split('|')[5]) + 2
                trans_gene_Dict = {':'.join([fileds[0],str(i)]):'' for i in range(1,length + 1)}
                
    # 2.add to full positions for interactome
    for full_key in inter_gene_Dict.keys():
        if full_key not in inter_Dict:
            inter_Dict[full_key] = 0.0
        else:
            pass
    # 2.1add to full positions for interactome
    for full_key in trans_gene_Dict.keys():
        if full_key not in trans_Dict:
            trans_Dict[full_key] = 0.0
        else:
            pass

    # 3.merge two dict information
    for key,val in trans_Dict.items():
        if key in inter_Dict:
            combined_Dict[key] = [float(val),float(inter_Dict[key])]
        else:
            combined_Dict[key] = [float(val),0]
    
    # 4.sort dict according to positions
    combined_Dict_newkey = {int(key.split(':')[1]):val for key,val in combined_Dict.items()}
    combined_Dict_sorted = {i:combined_Dict_newkey[i] for i in sorted(combined_Dict_newkey)}

    # 5.smooth density ratio
    slidedDict = SlideWindowEnrich(combined_Dict_sorted)
    # 5.1 get prefix id
    id = [key for key in combined_Dict.keys()][0].split(':')[0]
    cdsStart = int(id.split('|')[3]);cdsEnd = int(id.split('|')[4])
    cdsLength = cdsEnd - cdsStart + 1
    # 5.2 add id to key
    tmp_Dict = {':'.join([id,str(key)]):val for key,val in slidedDict.items()}
    
    # 6.filter CDS region and assign to new coordinates
    final_Dict = {}
    for key,val in tmp_Dict.items():
        id,pos = key.split(':')
        if cdsStart <= int(pos) <= cdsEnd:
            new_pos = int(pos) - cdsStart + 1
            new_key = ':'.join([id,str(new_pos)])
            final_Dict[new_key] = val
        else:
            pass

    return(final_Dict)

# batch enrichment function
def batchEnrichment(inputIPFile,inputInputFile,geneList,ouputFile):
    geneList = geneList

    # 1.batch for gene
    ratioFullDict = {}
    for g in geneList:
        ratioDict = SingleGeneEnrichment(inputIPFile=inputIPFile,inputInputFile=inputInputFile,geneName=g)
        ratioFullDict.update(ratioDict)

    # 2.outputfile
    outFile = open(ouputFile,'w')
    for key,val in ratioFullDict.items():
        id,pos = key.split(':')
        outFile.write('\t'.join([id,pos,str(val)]) + '\n')

    outFile.close()