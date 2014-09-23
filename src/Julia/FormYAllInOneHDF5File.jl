

using ArgParse
using HDF5

#Parse arguments function
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--out_file", "-o"
			help = "Output text file, on same basis as --file_names"
        "--file_names", "-f"
			help = "Text file of training file name prefixes: '\$(file_names[i])-\$(kmer_size)mers.kcount'"
        "--counts_HDF5_file", "-c"
			help = "Location of Kcounts in one big HDF5 file, datasets in the form '\$(file_names[i])-\$(kmer_size)mers.kcount'"
		"--input_file", "-i"
			help = "File name of sample kcount file"
        "--kmer_size", "-k"
			help = "Kmer size"
    end
    return parse_args(s)
end

#countPair function
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
	#Read in the args
	@everywhere parsed_args = remotecall_fetch(1,()->parse_commandline())
	println("Parsed args:")
	for (arg,val) in parsed_args
	    println("  $arg  =>  $val")
	end
	if parsed_args["out_file"]==nothing
	    error("Output file is required, use --help to see usage")
	end
	if parsed_args["file_names"]==nothing
	    error("Text file of file name prefixes is required, use --help to see usage")
	end
	if parsed_args["counts_HDF5_file"]==nothing
	    error("Must provide the counts file, use --help to see usage")
	end
	if parsed_args["kmer_size"]==nothing
	    error("kmer_size must be provided, use --help to see usage")
	end
	if parsed_args["input_file"]==nothing
	    error("input file must be provided, use --help to see usage")
	end
	#Use the command line arguments to populate variables
	@everywhere out_file = parsed_args["out_file"]
	@everywhere kmer_size = int(parsed_args["kmer_size"])
	@everywhere fid = open(parsed_args["file_names"],"r")
	@everywhere file_names = map(x->strip(x),readlines(fid))
	@everywhere close(fid)
	@everywhere input_file = parsed_args["input_file"]
	@everywhere sample_kmers = HDF5.h5read(input_file,"/kmers");
	@everywhere to_divide = float(sum(sample_kmers[2,:]));
	@everywhere counts_HDF5_file = parsed_args["counts_HDF5_file"]
	@everywhere all_database_kmers = HDF5.h5open("$(counts_HDF5_file)","r");
	
	#Sample Y vector
	Y = SharedArray(Float64, (length(file_names),1), init=0);
	
	#Count them in parallel
	@sync begin
		@parallel for i=1:length(file_names)
			database_kmers = all_database_kmers["$(file_names[i])"][:,:];
			temp = countPair(database_kmers,sample_kmers);
			Y[i] = float(temp[2])./to_divide;
		end
	end
	
	#Convert the array
	Y2=convert(Array,Y)
	
	#Delete output if it's there
	if isfile(out_file)
		rm(out_file)
	end
	
	#Write output
	#h5write(out_file,"/common_kmers",Y2)
	fid = open(out_file, "w")
	for i=1:length(Y2)
		write(fid, "$(Y2[i])\n")
	end
	
end

main()