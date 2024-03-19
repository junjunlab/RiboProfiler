using DataStructures

##############################################################################################################
# for mapping to genome qc analysis
##############################################################################################################

# define function
function MetageneAnalysis(;geneInfo,inputFile,outputFile,mode="st",type="codon",cdslength=600,expression=75,exclude=90)
    ################################################################
    # 1.choose which mode to analysis
    ################################################################
    # choose which mode to analysis
    if mode == "st"
        mylist = range(-50, 1500)
        # codon coorodinate
        region = [i for i in -51:3:1500]
        codonPos = [i for i in range(-17,500)]
    elseif mode == "sp"
        mylist = range(-1500, 50)
        # codon coorodinate
        region = [i for i in -1501:3:50]
        codonPos = [i for i in range(-500,17)]
    else
        println("pelase give the st/sp mode")
    end

    # save info
    rangeDict = Dict([i, 0.0] for i in mylist)
    countDict = Dict([i, 0.0] for i in mylist)

    ################################################################
    # 2.extract gene feature info
    ################################################################
    featureLenDict = Dict{String,Tuple}()

    # open gene info file
    open(geneInfo,"r") do input
        for line in eachline(input)
            fileds = split(line)

            # split tags
            _,_,gene_id,trans_id,_,_,_,_,utr5,cds,utr3 = split(line)
            # exon length
            utr5,cds,utr3 = parse(Int64,utr5),parse(Int64,cds),parse(Int64,utr3)
            # cds start and end position
            cdsStart,cdsEnd = (utr5 + 1),(utr5 + cds - 2)

            # save in dict
            featureLenDict[trans_id] = (cdsStart,cdsEnd,cds)
        end
    end

    ################################################################
    # 3.filter gene CDS length and expression higher than threshold
    ################################################################
    gene_infoDict = Dict{String,Float64}()
    filtedGeneDict = Dict{String,Float64}()

    # open gene info file
    open(inputFile,"r") do input
        for line in eachline(input)
            fileds = split(line)

            # split tags
            trans_id = fileds[2]
            transpos = parse(Int64,fileds[3])
            rpm = parse(Float64,fileds[4])

            # filter CDS > 600 nt gene
            if haskey(featureLenDict,trans_id)
                cdsStart,cdsEnd,cdsLength = featureLenDict[trans_id]
                if cdsLength >= cdslength && cdsLength%3 == 0
                    # not include first 30 codons to exclude high densitys around start codon
                    key = join([trans_id,cdsLength - exclude],":") # key
                    if !haskey(gene_infoDict,key)
                        gene_infoDict[key] = 0
                    else
                        if mode == "st"
                            if cdsStart + exclude <= transpos <= cdsEnd
                                gene_infoDict[key] += rpm
                            else
                                continue
                            end
                        else
                            if cdsStart <= transpos <= cdsEnd - exclude
                                gene_infoDict[key] += rpm
                            else
                                continue
                            end
                        end
                    end 
                else
                    continue
                end
            end
        end
    end

    # filter CDS expression >= 75
    for (key,val) in gene_infoDict
        if val >= expression
            meanNorm = val / parse(Int64,split(key,":")[2])
            filtedGeneDict[key] = meanNorm
        else
            continue
        end
    end

    ################################################################
    # 4.Meta-gene analysis from start codon
    ################################################################
    # open file
    open(inputFile,"r") do input
        for line in eachline(input)
            fileds = split(line)

            # split tags
            trans_id = fileds[2]
            transpos = parse(Int64,fileds[3])
            density = parse(Float64,fileds[4])
            cdsStart,cdsEnd,cdsLength = featureLenDict[trans_id]

            # key
            id = join([trans_id,cdsLength - exclude],":")

            # calculate -50-1500 sum denisty
            if haskey(filtedGeneDict,id)
                # calculate relative to start/stop codon distance
                if mode == "st"
                    reldist = transpos - cdsStart
                elseif mode == "sp"
                    reldist = transpos - cdsEnd
                else
                    println("pelase give the st/sp mode")
                end
                # sum up reads
                if mylist[1] <= reldist <= mylist[end]
                    # divide reads at one position by average number of reads per base for this gene
                    reads = density / filtedGeneDict[id]
                    rangeDict[reldist] += reads
                    # how often was this position counted in the calculation
                    countDict[reldist] += 1
                else
                    continue
                end
            else
                continue
            end
        end
    end

    ################################################################
    # 5.output data
    ################################################################
    fullDict = Dict{Int64,Float64}()
    for elem in zip(sort(rangeDict),sort(countDict))
        col0 = elem[1][1]       # list0 col0 = position (K)
        col1 = elem[1][2]       # list0 col1 = norm read number 
        col2 = elem[2][2]       # list1 col1 = how often was position counted

        #normalization2 by frequnecy
        if col2 == 0
            fullDict[col0] = 0
        else
            fullDict[col0] = col1 / col2
        end
    end

    # calculate relative density(position denisty/mean_density)

    new_fullDict = Dict{Float64,Float64}()
    meanDensity = sum(values(fullDict))/1551
    for (key,val) in fullDict
        relDensity = val/meanDensity
        new_fullDict[key] = relDensity
    end

    # condon position transform
    if type == "codon"
        if mode == "st"
            codonDict = Dict{Int64,Float64}()
            for i in range(1,length(region))
                codonRegion = range(region[i],region[i] + 2)
                count = 0
                codonDict[codonPos[i]] = 0
                # sum up codon three position densitys
                for j in codonRegion
                    if j in keys(new_fullDict)
                        codonDict[codonPos[i]] += new_fullDict[j]
                        count += 1
                    end
                end
                # codon mean density
                codonDict[codonPos[i]] = codonDict[codonPos[i]]/count
            end
            finalDict = codonDict
        elseif  mode == "sp"
            codonDict = Dict{Int64,Float64}()
            for i in range(1,length(region))
                codonRegion = range(region[i],region[i] + 2)
                count = 0
                codonDict[codonPos[i]] = 0
                # sum up codon three position densitys
                for j in codonRegion
                    if j in keys(new_fullDict)
                        codonDict[codonPos[i]] += new_fullDict[j]
                        count += 1
                    end
                end
                # codon mean density
                codonDict[codonPos[i]] = codonDict[codonPos[i]]/count
            end
            finalDict = codonDict
        else
            finalDict = new_fullDict
        end
    # nt position mode
    elseif type == "nt"
        finalDict = new_fullDict
    else
        print("pelase give the codon/nt mode")
    end
                
    # Finish output
    tupledlist = sort(finalDict)

    # output
    outFileP = open(outputFile, "w")
        
    for (pos,meanDensity) in tupledlist
        write(outFileP,join([pos,meanDensity],"\t")*"\n")
    end
    close(outFileP)
end


##############################################################################################################
# for mapping to transcriptome qc analysis
##############################################################################################################

# define function
function MetageneAnalysis_ontrans(;inputFile,outputFile,mode="st",type="codon",cdslength=600,expression=75,exclude=90)
    # choose which mode to analysis
    if mode == "st"
        mylist = range(-50, 1500)
        # codon coorodinate
        region = [i for i in -51:3:1500]
        codonPos = [i for i in range(-17,500)]
    elseif mode == "sp"
        mylist = range(-1500, 50)
        # codon coorodinate
        region = [i for i in -1501:3:50]
        codonPos = [i for i in range(-500,17)]
    else
        println("pelase give the st/sp mode")
    end

    # save info
    rangeDict = Dict([i, 0.0] for i in mylist)
    countDict = Dict([i, 0.0] for i in mylist)

    #########################################################
    # filter gene CDS length and expression higher thean threshold
    gene_infoDict = Dict{String,Float64}()
    filtedGeneDict = Dict{String,Float64}()
    
    # open file
    open(inputFile,"r") do input
        for line in eachline(input)
            fileds = split(line)

            # tags
            gene_name,_,trans_id,cdsStart,cdsEnd,_ = split(fileds[1],"|")
            pos = parse(Int64,fileds[2])
            cdsLength = parse(Int64,cdsEnd) - parse(Int64,cdsStart)
            density = parse(Float64,fileds[4])

            # filter CDS > 400 nt gene
            if cdsLength > cdslength && cdsLength%3 == 0
                # not include first 30 codons to exclude high densitys around start codon
                key = join([trans_id,cdsLength - exclude],":") # key
                if !haskey(gene_infoDict,key)
                    gene_infoDict[key] = 0
                else
                    if mode == "st"
                        if parse(Int64,cdsStart) + exclude <= pos <= parse(Int64,cdsEnd)
                            gene_infoDict[key] += density
                        else
                            continue
                        end
                    else
                        if parse(Int64,cdsStart) <= pos <= parse(Int64,cdsEnd) - exclude
                            gene_infoDict[key] += density
                        else
                            continue
                        end
                    end
                end
                
            else
                continue
            end
        end
    end

    # filter CDS expression > 50
    for (key,val) in gene_infoDict
        if val > expression
            meanNorm = val / parse(Int64,split(key,":")[2])
            filtedGeneDict[key] = meanNorm
        else
            continue
        end

    end
    #########################################################
    # Meta-gene analysis from start codon

    # open file
    open(inputFile,"r") do input
        for line in eachline(input)
            fileds = split(line)

            # tags
            gene_name,_,trans_id,cdsStart,cdsEnd,_ = split(fileds[1],"|")
            pos = parse(Int64,fileds[2])
            cdsLength = parse(Int64,cdsEnd) - parse(Int64,cdsStart)
            density = parse(Float64,fileds[4])
            id = join([trans_id,cdsLength - 90],":")

            # calculate -50-1500 sum denisty
            if haskey(filtedGeneDict,id)
                # calculate relative to start/stop codon distance
                if mode == "st"
                    reldist = pos - parse(Int64,cdsStart)
                elseif mode == "sp"
                    reldist = pos - parse(Int64,cdsEnd)
                else
                    println("pelase give the st/sp mode")
                end
                # sum up reads
                if mylist[1] <= reldist <= mylist[end]
                    # divide reads at one position by average number of reads per base for this gene
                    reads = density / filtedGeneDict[id]
                    rangeDict[reldist] += reads
                    # how often was this position counted in the calculation
                    countDict[reldist] += 1
                else
                    continue
                end
            else
                continue
            end
        end
    end
    
    #########################################################
    # output data

    fullDict = Dict{Int64,Float64}()
    for elem in zip(sort(rangeDict),sort(countDict))
        col0 = elem[1][1]      # list0 col0 = position (K)
        col1 = elem[1][2]       # list0 col1 = norm read number 
        col2 = elem[2][2]       # list1 col1 = how often was position counted

        #normalization2 by frequnecy
        if col2 == 0
            fullDict[col0] = 0
        else
            fullDict[col0] = col1 / col2
        end
    end

    # calculate relative density(position denisty/mean_density)

    new_fullDict = Dict{Float64,Float64}()
    meanDensity = sum(values(fullDict))/1551
    for (key,val) in fullDict
        relDensity = val/meanDensity
        new_fullDict[key] = relDensity
    end

    # condon position transform
    if type == "codon"
        if mode == "st"
            codonDict = Dict{Int64,Float64}()
            for i in range(1,length(region))
                codonRegion = range(region[i],region[i] + 2)
                count = 0
                codonDict[codonPos[i]] = 0
                # sum up codon three position densitys
                for j in codonRegion
                    if j in keys(new_fullDict)
                        codonDict[codonPos[i]] += new_fullDict[j]
                        count += 1
                    end
                end
                # codon mean density
                codonDict[codonPos[i]] = codonDict[codonPos[i]]/count
            end
            finalDict = codonDict
        elseif  mode == "sp"
            codonDict = Dict{Int64,Float64}()
            for i in range(1,length(region))
                codonRegion = range(region[i],region[i] + 2)
                count = 0
                codonDict[codonPos[i]] = 0
                # sum up codon three position densitys
                for j in codonRegion
                    if j in keys(new_fullDict)
                        codonDict[codonPos[i]] += new_fullDict[j]
                        count += 1
                    end
                end
                # codon mean density
                codonDict[codonPos[i]] = codonDict[codonPos[i]]/count
            end
            finalDict = codonDict
        else
            finalDict = new_fullDict
        end
    # nt position mode
    elseif type == "nt"
        finalDict = new_fullDict
    else
        print("pelase give the codon/nt mode")
    end
                
    # Finish output
    tupledlist = sort(finalDict)

    # output
    outFileP = open(outputFile, "w")
        
    for (pos,meanDensity) in tupledlist
        write(outFileP,join([pos,meanDensity],"\t")*"\n")
    end
    close(outFileP)
end