# ==============================================================================
# FormTrainingMatrix.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Forms the training matrices required for the CommonKmers method.
# ==============================================================================
# call with something like: julia -p 2 FormTrainingMatrix.jl -i ../Data/JFfiles.txt -o ../Output/CommonKmersMatrix.h5 -c ./count_in_file -s 2
using HDF5
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
		"--input_files", "-i"
			help = "File containing full paths to the 30mer or 50mer jellyfish files."
		"--output_file", "-o"
			help = "output file name. Example: /path/to/CommonKmersData/CommonKmersMatrix-30mers.h5"
		"--count_in_file_binary", "-c"
			help = "Location of the count_in_file_binary. eg ~/bin/./count_in_file"
        "--size_of_chunk", "-s"
            help = "Matrix is computed in chunks, this determines the chunk size"
            arg_type = Int
            default = 500
    end
    return parse_args(s)
end

#Parse the args
@everywhere parsed_args = remotecall_fetch(1,()->parse_commandline())
@everywhere input_files = parsed_args["input_files"]
@everywhere output_file = parsed_args["output_file"]
@everywhere count_in_file_binary = parsed_args["count_in_file_binary"]
@everywhere chunk_size = parsed_args["size_of_chunk"]

#Check for errors
if parsed_args["input_files"]==nothing
    error("A file containing the full paths to the training genomes is required, use --help to see usage")
end
if parsed_args["output_file"]==nothing
    error("An output file is required, use --help to see usage")
end
if parsed_args["count_in_file_binary"]==nothing
    error("count_in_file_binary location is required, use --help to see usage")
end

#Read in the file names to be processed
@everywhere fid = open(input_files,"r")
@everywhere file_names = split(readall(fid))
close(fid)
@everywhere num_seqs = length(file_names)

#Form the CommonKmer matrix
CommonKmersMatrix = SharedArray(Int64, (num_seqs,num_seqs), init=0);
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
					ikmers = "$(file_names[ii])"
					for jj=j:j+chunk_size-1
						if jj<=num_seqs
							jkmers = "$(file_names[jj])"
							(CommonKmersMatrix[ii,jj],CommonKmersMatrix[jj,ii]) = int64(split(readall(`$(count_in_file_binary) $(ikmers) $(jkmers)`)));
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
CommonKmersMatrix2=convert(Array,CommonKmersMatrix)'

if isfile(output_file)
	rm(output_file)
end

#Write the output file
h5write(output_file,"/common_kmers",CommonKmersMatrix2)













