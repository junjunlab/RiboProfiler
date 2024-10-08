# load package
using XAM
# using Pkg

##############################################################################################################
# for mapping to genome qc analysis
##############################################################################################################

###################################################################
# 1.load gene fearures positions and transform into trans position
###################################################################
function prepareQCdata(;longestTransInfo,inFile,outFile,seqType,assignType)
    geneinfoDict = Dict{String, Array}()

    open(longestTransInfo, "r") do geneinfo
        for line in eachline(geneinfo)
            
            # split tags
            _, gene_name, _, transid, chr, strand, cdsRg, exon, utr5, cds, utr3 = split(line)
            # exon length
            utr5, cds, utr3 = parse(Int64, utr5), parse(Int64, cds), parse(Int64, utr3)
            exonLength = utr5 + cds + utr3
            # cds start and end position
            cdsStart, cdsEnd = (utr5 + 1), (utr5 + cds - 2)
            
            # get exon positions
            exonPositions = []
            for rg in split(exon, ",")
                posSt, posEnd = parse(Int64, split(rg, ":")[1]), parse(Int64, split(rg, ":")[2])
                push!(exonPositions, collect(posSt:posEnd))
            end
            exonPositions = reduce(vcat, exonPositions)
            
            # calculate geneinfoDict using vectorized operations
            counts = 1:exonLength
            if strand == "+"
                counts = counts'
            else
                counts = reverse(counts')
            end
            
            posKeys = ["$chr|$pos" for pos in exonPositions]

            for i in range(1,length(posKeys))
                geneinfoDict[posKeys[i]] = [cdsStart, cdsEnd,counts[i],transid]
            end
        end
    end
    
    display("Transforming genomic positions into transcriptome positions has been done successfully.")

    ################################################################
    # 2.calculate frame on each gene and distance to start/stop codon
    ################################################################

    # define function
    function RiboQcAnalysis(inputFile,outputFile;seq_type,assignType)
        
        # save in dict
        frame_dict = Dict{String,Int64}()
        # total_mapped_counts = 0

        # open sam file
        reader = open(SAM.Reader,inputFile)
        record = SAM.Record()

        # loop
        while !eof(reader)
            empty!(record)
            read!(reader, record)
        # for record in reader
            # do something
            if SAM.ismapped(record) # (remove flag4)
                # tags
                refname,align_pos,read_length = SAM.refname(record),SAM.position(record),SAM.seqlength(record)

                # read flag tag
                flag = BAM.flag(record)

                # flag16(+) use 5'end as alignpos and flag0(-) use 3'end as alignpos
                if seq_type == "singleEnd"
                    # flag 0 reads from - stand
                    if flag == 16
                        if assignType == "end5"
                            # read key
                            pos = align_pos + read_length - 1
                            readKey = join([refname,pos],"|")
                        else
                            readKey = join([refname,align_pos],"|")
                        end
                    elseif flag == 0
                        # flag 16 reads from + stand
                        if assignType == "end5"
                            # read key
                            readKey = join([refname,align_pos],"|")
                        else
                            pos = align_pos + read_length - 1
                            readKey = join([refname,pos],"|")
                        end
                    else
                        println("There are other flags!")
                    end
                elseif seq_type == "pairedEnd"
                    if flag == 16
                        # flag 16 reads from + stand
                        if assignType == "end5"
                            # read key
                            readKey = join([refname,align_pos],"|")
                        else
                            pos = align_pos + read_length - 1
                            readKey = join([refname,pos],"|")
                        end
                    elseif flag == 0
                        # flag 0 reads from - stand
                        if assignType == "end5"
                            # read key
                            pos = align_pos + read_length - 1
                            readKey = join([refname,pos],"|")
                        else
                            readKey = join([refname,align_pos],"|")
                        end
                    else
                        println("There are other flags!")
                    end
                end

                # get relative distance from start/stop codon and frame information
                if haskey(geneinfoDict,readKey)

                    # get gene info
                    start_codon_pos,stop_codon_pos,transPos,transid = geneinfoDict[readKey]
                    
                    rel2st = transPos - start_codon_pos
                    rel2sp = transPos - stop_codon_pos

                    # assign frame
                    frame_st = abs(rel2st)%3
                    frame_sp = abs(rel2sp)%3
                    
                    # read center position
                    if assignType == "end5"
                        align_pos_center = transPos + (read_length ÷ 2)
                    elseif assignType == "end3"
                        align_pos_center = transPos - (read_length ÷ 2)
                    else
                        println("There are other flags!")
                    end

                    # feaure type(5UTR,CDS,3UTR)
                    if align_pos_center <= start_codon_pos
                        ftype = 2 # 5utr
                    elseif start_codon_pos < align_pos_center <= stop_codon_pos + 2
                        ftype = 3 # cds
                    else stop_codon_pos + 2 < align_pos_center
                        ftype = 1 # 3UTR
                    end

                    # key
                    key = join([read_length,frame_st,rel2st,frame_sp,rel2sp,ftype,transPos,transid],"\t")

                    # init dict and count
                    if !haskey(frame_dict,key)
                        frame_dict[key] = 1
                    else
                        frame_dict[key] += 1
                    end
                end
            end
        end

        close(reader)
        ################################################################
        # 3.output results
        ################################################################
        # output file
        outfile = open(outputFile,"w")
        for (key,val) in frame_dict
            # write(outfile,"$key\t$val\n")
            # normalized_counts = (val/total_mapped_counts)*1000000
            write(outfile,"$key\t$val\n")
        end
        close(outfile)
    end

    # loop for output
    display("Processing mapping files...")
    inFile = [s for s in split(inFile,",")]
    outFile = [o for o in split(outFile,",")]
    for i in range(1,length(inFile))
        RiboQcAnalysis(inFile[i],outFile[i],seq_type = seqType,assignType = assignType)
        tmp_name = inFile[i]
        display("$tmp_name has been processed.")
    end
end


###############################################################################################################
# XAM verson for 0.4.0
###############################################################################################################
function prepareQCdata2(;longestTransInfo,inFile,outFile,seqType,assignType)
    geneinfoDict = Dict{String, Array}()

    open(longestTransInfo, "r") do geneinfo
        for line in eachline(geneinfo)
            
            # split tags
            _, gene_name, _, transid, chr, strand, cdsRg, exon, utr5, cds, utr3 = split(line)
            # exon length
            utr5, cds, utr3 = parse(Int64, utr5), parse(Int64, cds), parse(Int64, utr3)
            exonLength = utr5 + cds + utr3
            # cds start and end position
            cdsStart, cdsEnd = (utr5 + 1), (utr5 + cds - 2)
            
            # get exon positions
            exonPositions = []
            for rg in split(exon, ",")
                posSt, posEnd = parse(Int64, split(rg, ":")[1]), parse(Int64, split(rg, ":")[2])
                push!(exonPositions, collect(posSt:posEnd))
            end
            exonPositions = reduce(vcat, exonPositions)
            
            # calculate geneinfoDict using vectorized operations
            counts = 1:exonLength
            if strand == "+"
                counts = counts'
            else
                counts = reverse(counts')
            end
            
            posKeys = ["$chr|$pos" for pos in exonPositions]

            for i in range(1,length(posKeys))
                geneinfoDict[posKeys[i]] = [cdsStart, cdsEnd,counts[i],transid]
            end
        end
    end
    
    display("Transforming genomic positions into transcriptome positions has been done successfully.")

    ################################################################
    # 2.calculate frame on each gene and distance to start/stop codon
    ################################################################

    # define function
    function RiboQcAnalysis(inputFile,outputFile;seq_type,assignType)
        
        # save in dict
        frame_dict = Dict{String,Int64}()
        # total_mapped_counts = 0

        # open sam file
        reader = open(SAM.Reader,inputFile)
        record = SAM.Record()

        # loop
        while !eof(reader)
            empty!(record)
            read!(reader, record)
        # for record in reader
            # do something
            if SAM.ismapped(record) # (remove flag4)
                # tags
                refname,align_pos,read_length = SAM.refname(record),SAM.position(record),SAM.seqlength(record)

                # read flag tag
                flag = BAM.flags(record)

                # flag16(+) use 5'end as alignpos and flag0(-) use 3'end as alignpos
                if seq_type == "singleEnd"
                    # flag 0 reads from - stand
                    if flag == 16
                        if assignType == "end5"
                            # read key
                            pos = align_pos + read_length - 1
                            readKey = join([refname,pos],"|")
                        else
                            readKey = join([refname,align_pos],"|")
                        end
                    elseif flag == 0
                        # flag 16 reads from + stand
                        if assignType == "end5"
                            # read key
                            readKey = join([refname,align_pos],"|")
                        else
                            pos = align_pos + read_length - 1
                            readKey = join([refname,pos],"|")
                        end
                    else
                        println("There are other flags!")
                    end
                elseif seq_type == "pairedEnd"
                    if flag == 16
                        # flag 16 reads from + stand
                        if assignType == "end5"
                            # read key
                            readKey = join([refname,align_pos],"|")
                        else
                            pos = align_pos + read_length - 1
                            readKey = join([refname,pos],"|")
                        end
                    elseif flag == 0
                        # flag 0 reads from - stand
                        if assignType == "end5"
                            # read key
                            pos = align_pos + read_length - 1
                            readKey = join([refname,pos],"|")
                        else
                            readKey = join([refname,align_pos],"|")
                        end
                    else
                        println("There are other flags!")
                    end
                end

                # get relative distance from start/stop codon and frame information
                if haskey(geneinfoDict,readKey)

                    # get gene info
                    start_codon_pos,stop_codon_pos,transPos,transid = geneinfoDict[readKey]
                    
                    rel2st = transPos - start_codon_pos
                    rel2sp = transPos - stop_codon_pos

                    # assign frame
                    frame_st = abs(rel2st)%3
                    frame_sp = abs(rel2sp)%3
                    
                    # read center position
                    if assignType == "end5"
                        align_pos_center = transPos + (read_length ÷ 2)
                    elseif assignType == "end3"
                        align_pos_center = transPos - (read_length ÷ 2)
                    else
                        println("There are other flags!")
                    end

                    # feaure type(5UTR,CDS,3UTR)
                    if align_pos_center <= start_codon_pos
                        ftype = 2 # 5utr
                    elseif start_codon_pos < align_pos_center <= stop_codon_pos + 2
                        ftype = 3 # cds
                    else stop_codon_pos + 2 < align_pos_center
                        ftype = 1 # 3UTR
                    end

                    # key
                    key = join([read_length,frame_st,rel2st,frame_sp,rel2sp,ftype,transPos,transid],"\t")

                    # init dict and count
                    if !haskey(frame_dict,key)
                        frame_dict[key] = 1
                    else
                        frame_dict[key] += 1
                    end
                end
            end
        end

        close(reader)
        ################################################################
        # 3.output results
        ################################################################
        # output file
        outfile = open(outputFile,"w")
        for (key,val) in frame_dict
            # write(outfile,"$key\t$val\n")
            # normalized_counts = (val/total_mapped_counts)*1000000
            write(outfile,"$key\t$val\n")
        end
        close(outfile)
    end

    # loop for output
    display("Processing mapping files...")
    inFile = [s for s in split(inFile,",")]
    outFile = [o for o in split(outFile,",")]
    for i in range(1,length(inFile))
        RiboQcAnalysis(inFile[i],outFile[i],seq_type = seqType,assignType = assignType)
        tmp_name = inFile[i]
        display("$tmp_name has been processed.")
    end
end