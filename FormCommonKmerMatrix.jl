#!/usr/local/julia/usr/bin/julia
#Example usage
#julia -p 2 test.jl --counts_dir /raid2/labs/Koslicki_lab/koslickd/CommonKmers/TrainingOnRepoPhlAn/Counts/Counts20C/ --size_of_chunk 2 --file_names testFileNames.txt  -o test_out.h5 -k 20

using ArgParse
using HDF5

#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--counts_dir", "-c"
            help = "Path to the directory that contains the *.kcount files"
        "--size_of_chunk", "-s"
            help = "Matrix is computed in chunks, this determines the chunk size"
            arg_type = Int
            default = 500
        "--file_names", "-f"
            help = "Text file of file name prefixes: '\$(file_names[i])-\$(kmer_size)mers.kcount'"
        "--output_file", "-o"
            help = "Full name of the output HDF5 file"
        "--kmer_size", "-k"
            help = "Kmer size"
    end
    return parse_args(s)
end

#Define the countsPair function
@everywhere function countPair(iKmers, jKmers)
    iMax = size(iKmers)[2]
    jMax = size(jKmers)[2]
    i = 1
    j = 1
    iTotal = 0
    jTotal = 0
    while (i < iMax) & (j < jMax)
        if iKmers[1,i] == jKmers[1,j]
            iTotal += iKmers[2,i]
            jTotal += jKmers[2,j]
            i += 1
            j += 1
        elseif iKmers[1,i] < jKmers[1,j]
            i += 1
        else
            j += 1
        end
    end
    return(iTotal, jTotal)
end

function main()
	#Read in the command line arguments and error check them
	@everywhere parsed_args = remotecall_fetch(1,()->parse_commandline())
	println("Parsed args:")
	for (arg,val) in parsed_args
	    println("  $arg  =>  $val")
	end
	if parsed_args["output_file"]==nothing
	    error("An output file is required, use --help to see usage")
	end
	if parsed_args["file_names"]==nothing
	    error("Text file of file name prefixes is required, use --help to see usage")
	end
	if parsed_args["counts_dir"]==nothing
	    error("Must provide the counts directory, use --help to see usage")
	end
	if parsed_args["kmer_size"]==nothing
	    error("kmer_size must be provided, use --help to see usage")
	end
	#Use the command line arguments to populate variables
	@everywhere counts_dir = parsed_args["counts_dir"]
	@everywhere kmer_size = int(parsed_args["kmer_size"])
	@everywhere fid = open(parsed_args["file_names"],"r")
	@everywhere file_names = split(readall(fid))
	close(fid)
	@everywhere num_seqs = length(file_names)
	@everywhere out_file = parsed_args["output_file"]
	@everywhere chunk_size = int(parsed_args["size_of_chunk"])


	#This is the common kmers matrix
	CommonKmersMatrix = SharedArray(Int64, (num_seqs,num_seqs), init=0);

	#Populate the matrix in parallel, in chunks
	for j = 1:chunk_size:num_seqs
		for i = 1:chunk_size:num_seqs
			#only do the subdiagonal
			if i>=j
				tic()
				print("On (row,column) chunk: ($(i),$(j)) of $(num_seqs)\n")
				#parallelize over the chunk
				@sync begin
				@parallel for ii = i:i+chunk_size-1
					if ii<=num_seqs
						iKmers = HDF5.h5read(string(counts_dir,"$(file_names[ii])-$(kmer_size)mers.kcount"),"/kmers");
						for jj=j:j+chunk_size-1
							if jj<=num_seqs
								jKmers = HDF5.h5read(string(counts_dir,"$(file_names[jj])-$(kmer_size)mers.kcount"),"/kmers");
								(CommonKmersMatrix[ii,jj],CommonKmersMatrix[jj,ii]) = countPair(iKmers,jKmers);
							end
						end
					end
				end
				end
				toc()
				gc()
			end
		end
	end

	#Convert the matrix into something that HDF5 can write
	CommonKmersMatrix2=convert(Array,CommonKmersMatrix)

	if isfile(out_file)
		rm(out_file)
	end

	h5write(out_file,"/common_kmers",CommonKmersMatrix2)
end

main()