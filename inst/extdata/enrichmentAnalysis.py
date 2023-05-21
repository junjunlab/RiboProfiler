# define function
def enrichmentAnalysis(IPfile,Inputfile,outputFile):
    inter_Dict = {}
    trans_Dict = {}
    combined_Dict = {}

    # 1.open denisty file and save in dict
    with open(IPfile,'r') as interactome, \
        open(Inputfile,'r') as translatome:
        # opne IP file
        for line in interactome:
            fileds = line.split()
            key = ':'.join([fileds[0],str(fileds[1])])
            inter_Dict[key] = fileds[3]
        # opne Input file
        for line in translatome:
            fileds = line.split()
            key = ':'.join([fileds[0],str(fileds[1])])
            trans_Dict[key] = fileds[3]

    # 2.merge two dict information
    for key,val in trans_Dict.items():
        if key in inter_Dict:
            combined_Dict[key] = [float(val),float(inter_Dict[key])]
        else:
            combined_Dict[key] = [float(val),0]

    # 3.get ratio
    for key,val in combined_Dict.items():
        if val[0] == 0:
            ratio = 0
            combined_Dict[key].append(ratio)
        else:
            ratio = val[1]/val[0]
            combined_Dict[key].append(ratio)

    # 4.output
    outFileP = open(outputFile, 'w')
        
    for key,val in combined_Dict.items():
        elem = key.split(':')
        outFileP.write('\t'.join([elem[0],str(elem[1]),str(val[0]),str(val[1]),str(val[2])]) + '\n')
    outFileP.close()