#!/usr/local/julia/usr/bin/julia

using ArgParse
using HDF5

#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--counts_dir", "-c"
		help = "Path to the directory that will contain the *.kcount files"
        "--file_names", "-f"
		help = "Text file of file name prefixes: '\$(file_names[i])-\$(kmer_size)mers.kcount'"
        "--output_file", "-o"
            help = "Output file"
        "--level_of_compression", "-l"
            help = "Compression level 0-9"
        "--kmer_size", "-k"
            help = "kmer size"
    end
    return parse_args(s)
end

function main()
	#Read in the args
	parsed_args = parse_commandline()
	println("Parsed args:")
	for (arg,val) in parsed_args
	    println("  $arg  =>  $val")
	end
	if parsed_args["counts_dir"]==nothing
	    error("Output directory is required, use --help to see usage")
	end
	if parsed_args["file_names"]==nothing
	    error("Text file of file name prefixes is required, use --help to see usage")
	end
	if parsed_args["output_file"]==nothing
	    error("Must provide the sequences directory, use --help to see usage")
	end
	if parsed_args["level_of_compression"]==nothing
	    error("Compression level must be provided, use --help to see usage")
	end
	if parsed_args["kmer_size"]==nothing
	    error("kmer size must be provided, use --help to see usage")
	end

	#populate the args
	counts_dir = parsed_args["counts_dir"]
	output_file = parsed_args["output_file"]
	level_of_compression = int(parsed_args["level_of_compression"])
	kmer_size = int(parsed_args["kmer_size"])
	fid = open(parsed_args["file_names"],"r")
	file_names = split(readall(fid))
	close(fid)
	num_seqs = length(file_names)

	if isfile(output_file)
		rm(output_file)
	end
	
	AllKmers_fid = h5open(output_file,"w")

	for i=1:num_seqs
		kmers=h5read(string(counts_dir,file_names[i],"-$(kmer_size)mers.kcount"),"/kmers")
		if level_of_compression == 0
			write(AllKmers_fid, convert(ASCIIString, file_names[i]), kmers)
		else
			AllKmers_fid[convert(ASCIIString, file_names[i]), "chunk", size(kmers), "compress", level_of_compression] = kmers
		end
	end

	close(AllKmers_fid)
end

main()