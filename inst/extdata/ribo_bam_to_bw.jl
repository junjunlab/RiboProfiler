using XAM
using BigWig

function ribo_bam_to_bw(;bam_file,output_file,seq_type,assignType,normalization,min_length,max_length,offset)
    ##############################################################################################################
    # read bam header and fetch chromsome size
    ##############################################################################################################
    reader = open(BAM.Reader,bam_file)
    
    chromsome_size = Dict()
    header = findall(BAM.header(reader), "SQ")

    # save chromosome info in dict
    for h in header
        tmp = split(string(h),"\t")
        chr = split(tmp[2],":")[2]
        len = parse(Int,split(tmp[3],":")[2])

        # save
        chromsome_size[chr] = len
    end

    ##############################################################################################################
    # process read offset
    ##############################################################################################################
    if length(offset) > 1
        offset_dict = Dict()

        fileds = split(offset,"|")
        read_len_select = [parse(Int,i) for i in split(fileds[1],",")] 
        read_offset_select = [parse(Int,i) for i in split(fileds[2],",")]

        # check read and offset length
        if length(read_len_select) == length(read_offset_select)
            for i in range(1,4)
                offset_dict[read_len_select[i]] = read_offset_select[i]
            end
        else
            println("Please make sure reads length is equal to offset!")
        end
        
    end

    ##############################################################################################################
    # process bam file
    ##############################################################################################################
    read_dict = Dict()

    # open sam file
    reader = open(BAM.Reader,bam_file)
    record = BAM.Record()

    # loop
    while !eof(reader)
        empty!(record)
        read!(reader, record)
        # do something
        if BAM.ismapped(record)
            # tags
            refname,align_pos,read_length = BAM.refname(record),BAM.position(record),BAM.seqlength(record)

            # check read length
            if min_length <= read_length <= max_length
                # read flag tag
                flag = BAM.flag(record)
                
                # flag16(+) use 5'end as alignpos and flag0(-) use 3'end as alignpos
                if seq_type == "singleEnd"
                    if flag == 16
                        # flag 0 reads from - stand
                        if assignType == "end5"
                            pos = align_pos + read_length - 1
                        else
                            pos = align_pos
                        end
                        # shift read length
                        if length(offset) > 0
                            pos = pos - offset_dict[read_length]
                        end
                    elseif flag == 0
                        # flag 16 reads from + stand
                        if assignType == "end5"
                            pos = align_pos
                        else
                            pos = align_pos + read_length - 1
                        end
                        # shift read length
                        if length(offset) > 0
                            pos = pos + offset_dict[read_length]
                        end
                    else
                        println("There are other flags!")
                    end
                elseif seq_type == "pairedEnd"
                    if flag == 16
                        # flag 16 reads from + stand
                        if assignType == "end5"
                            pos = align_pos
                        else
                            pos = align_pos + read_length - 1
                        end
                        # shift read length
                        if length(offset) > 0
                            pos = pos + offset_dict[read_length]
                        end
                    elseif flag == 0
                        # flag 0 reads from - stand
                        if assignType == "end5"
                            pos = align_pos + read_length - 1
                        else
                            pos = align_pos
                        end
                        # shift read length
                        if length(offset) > 0
                            pos = pos - offset_dict[read_length]
                        end
                    else
                        println("There are other flags!")
                    end
                end

                # save
                key = join([refname,pos],"|")

                if !haskey(read_dict,key)
                    read_dict[key] = 1
                else
                    read_dict[key] += 1
                end
            end
        end
    end

    ##############################################################################################################
    # normalization tpm or raw counts
    ##############################################################################################################

    if normalization == "rpm"
        total_counts = sum(values(read_dict))
        total_counts

        for (key,val) in read_dict
            norm_val = (val/total_counts)*1000000
            read_dict[key] = norm_val
        end
    end

    
    # sort chromosome and positions
    processed_keys = [(split(key, "|")[1], parse(Int, split(key, "|")[2]), key) for key in keys(read_dict)]
    sorted_entries = sort(processed_keys, by = x -> (x[1], x[2]))

    ##############################################################################################################
    # normalization tpm or raw counts
    ##############################################################################################################
    output_bw = open(output_file, "w")

    writer = BigWig.Writer(output_bw, [(chr, Integer(size)) for (chr,size) in chromsome_size])

    for (chr, pos, key) in sorted_entries
        counts = Float64(read_dict[key])
        write(writer, (String(chr), pos, pos, counts))
    end
    close(writer)
end