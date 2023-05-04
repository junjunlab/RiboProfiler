# load package
using XAM

# define function
function CalculateRibosomeDensity(;inputFile,outputFile,min = 23,max = 35)
    # save in dict
    density_dict = Dict{String,Float64}()

    # open sam file
    reader = open(SAM.Reader,inputFile)
    record = SAM.Record()

    # loop
    while !eof(reader)
        empty!(record)
        read!(reader, record)
        # do something
        if SAM.ismapped(record)
            # tags
            refname = SAM.refname(record)
            align_pos = SAM.position(record)
            read_length = SAM.seqlength(record)

            # filter read length
            if min <= read_length <= max
                end5 = align_pos
                end3 = end5 + read_length - 1
                # shift +- 11nt
                centerEnd5 = end5 + 11
                centerEnd3 = end3 - 11
                centerLength = centerEnd3 - centerEnd5 + 1
                
                # ribo density
                for elem in range(centerEnd5, centerEnd3)
                    key = "$refname:$elem"
                    if !haskey(density_dict,key)
                        density_dict[key] = (1.0 / centerLength)
                    else
                        density_dict[key] += (1.0 / centerLength)
                    end
                end
            end    
        end
    end

    # output file
    outfile = open(outputFile,"w")

    # total densitys
    total_density = sum(values(density_dict))

    # sort keys
    for key in sort(collect(keys(density_dict)))
        id,align_pos = split(key,":")
        raw_denisty = density_dict[key]
        
        # RPM normalization
        rpm = (raw_denisty/total_density)*1000000
        write(outfile,"$id\t$align_pos\t$raw_denisty\t$rpm\n")
    end
    close(outfile)
end