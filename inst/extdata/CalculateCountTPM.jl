# load package
using XAM

# define function
function CalculateCountTPM(;inputFile,outputFile,inputType)
    #####################################################
    # 1.count gene numbers
    #####################################################

    # save in dict
    gene_dict = Dict{String,Int64}()

    # open sam file
    reader = open(SAM.Reader,inputFile)
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

            # tags
            gene_name,_,_,cdsST,cdsSP,gene_length = split(refname,"|")
            
            if inputType == "ribo"
                # tags
                align_pos = SAM.position(record)
                read_length = SAM.seqlength(record)

                # read center position
                align_pos_center = align_pos + (read_length ÷ 2)
                
                # cds length
                cdsST,cdsSP = parse(Int64,cdsST),parse(Int64,cdsSP)
                cdsLength = cdsSP - cdsST + 1

                key = "$gene_name|$cdsLength"
                # count reads ribo on CDS region
                if cdsST <= align_pos_center <= cdsSP
                    if !haskey(gene_dict,key)
                        gene_dict[key] = 1
                    else
                        gene_dict[key] += 1
                    end
                end
            elseif inputType == "rna"
                # count reads rna on transcript region
                key = "$gene_name|$gene_length"
                if !haskey(gene_dict,key)
                    gene_dict[key] = 1
                else
                    gene_dict[key] += 1
                end
            else
                println("error!")
                break
            end
        end
    end
    close(reader)
    
    #####################################################
    # 2.calculate TPM values
    #####################################################

    # total reads one sample

    tpm_Dict = Dict{String,Array}()

    # get TPM values
    for (key,val) in gene_dict
        gene_name,gene_length = split(key,"|")

        # reads / geneLength(kb)
        rpk = val / (parse(Int64,gene_length) / 1000)
        tpm_Dict[gene_name] = [val,rpk]
    end

    # output
    tpm_output = open(outputFile,"w")
    totol_rpk = sum([rpk for (count,rpk) in values(tpm_Dict)])

    for (key,val) in tpm_Dict
        tpm = (val[2] / totol_rpk)*1000000
        # save
        write(tpm_output,join([key,val[1],tpm],"\t")*"\n")
    end
    close(tpm_output)
end