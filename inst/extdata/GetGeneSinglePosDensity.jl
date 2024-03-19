 # define function
function GetGeneSinglePosDensity(;geneInfo,inputFile,outputFile)
    ################################################################
    # 1.load density file
    ################################################################
    # save in dict
    densityDict = Dict{String,Float64}()
    chrDict = Dict{String,Int64}()

    open(inputFile,"r") do Density
        for line in eachline(Density)
            # split tags
            chr,alginPos,_,rpm = split(line)
            densityDict[join([chr,alginPos],"|")] = parse(Float64,rpm)
            if !haskey(chrDict,chr)
                chrDict[chr] = 0
            else
                continue
            end
        end
    end

    ################################################################
    # 2.transform into transcriptome coordinate
    ################################################################
    outGeneDensity = open(outputFile,"w")

    open(geneInfo,"r") do geneinfo
        for line in eachline(geneinfo)
            # split tags
            _,gene_name,gene_id,trans_id,chr,strand,_,exon,utr5,cds,utr3 = split(line)
            # exon length
            utr5,cds,utr3 = parse(Int64,utr5),parse(Int64,cds),parse(Int64,utr3)

            # exclude no 5UTR genes and can't be divided exactly by 3 genes
            # if utr5 != 0 && cds%3 == 0
            exonLength = utr5 + cds +utr3
            
            # cds start and end position
            cdsStart,cdsEnd = (utr5 + 1),(utr5 + cds - 2)

            # loop for every exon regions 
            if strand == "+" # + strand gene
                count = 0
                val = 1
            else # - strand gene
                count = exonLength + 1
                val = -1
            end

            # fill denisty
            for rg in split(exon,",")
                posSt,posEnd = parse(Int64,split(rg,":")[1]),parse(Int64,split(rg,":")[2])
                for chrPos in range(posSt,posEnd)
                    count += val
                    if chr in keys(chrDict)
                        density = get(densityDict,join([chr,chrPos],"|"),0)
                        write(outGeneDensity,join([gene_name,trans_id,count,density],"\t")*"\n")
                    else
                        continue
                    end
                end
            end 
        end
    end

    # close file
    close(outGeneDensity)
end