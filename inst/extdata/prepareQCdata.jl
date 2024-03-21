# load package
using XAM

##############################################################################################################
# for mapping to genome qc analysis
##############################################################################################################

###################################################################
# 1.load gene fearures positions and transform into trans position
###################################################################
function prepareQCdata(;longestTransInfo,samFile,outFile,seqType)
    geneinfoDict = Dict{String, Vector{Int}}()

    open(longestTransInfo, "r") do geneinfo
        for line in eachline(geneinfo)
            
            # split tags
            _, gene_name, _, _, chr, strand, cdsRg, exon, utr5, cds, utr3 = split(line)
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
                geneinfoDict[posKeys[i]] = [cdsStart, cdsEnd,counts[i]]
            end
        end
    end
    
    display("Transforming genomic positions into transcriptome positions has been done successfully.")

    ################################################################
    # 2.calculate frame on each gene and distance to start/stop codon
    ################################################################

    # define function
    function RiboQcAnalysis(inputFile,outputFile;seq_type)
        
        # save in dict
        frame_dict = Dict{String,Int64}()

        # open sam file
        reader = open(SAM.Reader,inputFile)
        # record = SAM.Record()

        # loop
        # while !eof(reader)
        #     empty!(record)
        #     read!(reader, record)
        for record in reader
            # do something
            if SAM.ismapped(record) # (remove flag4)
                # tags
                refname,align_pos,read_length = SAM.refname(record),SAM.position(record),SAM.seqlength(record)

                # read flag tag
                flag = SAM.flags(record)

                # flag16(+) use 5'end as alignpos and flag0(-) use 3'end as alignpos
                if seq_type == "singleEnd"
                    if flag == 16
                        # flag 0 reads from - stand
                        end3Pos = align_pos + read_length - 1
                        # read key
                        readKey = join([refname,end3Pos],"|")
                    elseif flag == 0
                        # flag 16 reads from + stand
                        # read key
                        readKey = join([refname,align_pos],"|")
                    else
                        println("There are other flags!")
                    end
                elseif seq_type == "pairedEnd"
                    if flag == 16
                        # flag 16 reads from + stand
                        # read key
                        readKey = join([refname,align_pos],"|")
                    elseif flag == 0
                        # flag 0 reads from - stand
                        end3Pos = align_pos + read_length - 1
                        # read key
                        readKey = join([refname,end3Pos],"|")
                    else
                        println("There are other flags!")
                    end
                end

                # get relative distance from start/stop codon and frame information
                if haskey(geneinfoDict,readKey)

                    # get gene info
                    start_codon_pos,stop_codon_pos,transPos = geneinfoDict[readKey]
                    
                    # relative distance
                    rel2st = transPos - start_codon_pos
                    rel2sp = transPos - stop_codon_pos

                    # assign frame
                    frame_st = abs(rel2st)%3
                    frame_sp = abs(rel2sp)%3
                    
                    # read center position
                    if flag == 16
                        align_pos_center = transPos + (read_length ÷ 2)
                    elseif flag == 0
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
                    key = join([read_length,frame_st,rel2st,frame_sp,rel2sp,ftype],"\t")

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
            write(outfile,"$key\t$val\n")
        end
        close(outfile)
    end

    # loop for output
    display("Processing sam files...")
    samFile = [s for s in split(samFile,",")]
    outFile = [o for o in split(outFile,",")]
    for i in range(1,length(samFile))
        RiboQcAnalysis(samFile[i],outFile[i],seq_type = seqType)
        tmp_name = samFile[i]
        display("$tmp_name has been processed.")
    end
end


##############################################################################################################
# for mapping to transcriptome qc analysis
##############################################################################################################

# define function
function prepareQCdata_ontrans(;samFile,outFile,seqType)
    # save in dict
    frame_dict = Dict{String,Int64}()

    # open sam file
    reader = open(SAM.Reader,samFile)
    # record = SAM.Record()

    # loop
    # while !eof(reader)
    #     empty!(record)
    #     read!(reader, record)
    for record in reader
        # do something
        if SAM.ismapped(record)
            # tags
            refname = SAM.refname(record)
            align_pos = SAM.position(record)
            read_length = SAM.seqlength(record)
            flag = SAM.flags(record)

            # flag0(-strand gene) and flag16(+strand gene) for read1
            if seqType == "singleEnd"
                if flag == 0
                    exact_pos = align_pos
                elseif flag == 16
                    exact_pos = align_pos + read_length - 1
                else
                    println("There are other flags!")
                end
            elseif seqType == "pairedEnd"
                if flag == 0
                    exact_pos = align_pos + read_length - 1
                elseif flag == 16
                    exact_pos = align_pos
                else
                    println("There are other flags!")
                end
            end
            

            # start and stop codon position
            start_codon_pos = parse(Int,split(refname,"|")[4])
            stop_codon_pos = parse(Int,split(refname,"|")[5])

            # relative distance
            rel2st = exact_pos - start_codon_pos
            rel2sp = exact_pos - (stop_codon_pos + 1)

            # assign frame
            frame_st = abs(rel2st)%3
            frame_sp = abs(rel2sp)%3

            # read center position
            if flag == 16
                align_pos_center = exact_pos + (read_length ÷ 2)
            elseif flag == 0
                align_pos_center = exact_pos - (read_length ÷ 2)
            else
                println("There are other flags!")
            end

            # assign 5'UTR/CDS/3'UTR features
            if align_pos_center < start_codon_pos
                ftype = 2 # 5'UTR
            elseif start_codon_pos <= align_pos_center <= stop_codon_pos
                ftype = 3 # CDS
            else align_pos_center > stop_codon_pos
                ftype = 1 # 3'UTR
            end

            # key
            key = join([read_length,frame_st,rel2st,frame_sp,rel2sp,ftype],"\t")
            
            # init dict and count
            if !haskey(frame_dict,key)
                frame_dict[key] = 1
            else
                frame_dict[key] += 1
            end
        end
    end

    close(reader)
    # output file
    outfile = open(outFile,"w")
    for (key,val) in frame_dict
        write(outfile,"$key\t$val\n")
    end
    close(outfile)
end