#!/usr/local/julia/usr/bin/julia
#Example usage 
#julia -p 3 FormCountsAndKCounts.jl -o Test/ -f Test/testFileNames.txt -s ../microbes/ -k 21 -j /home/pi/koslickd/jellyfish-2.1.1/bin/./jellyfish

using ArgParse
using HDF5

#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--output_dir", "-o"
		help = "Path to the directory that will contain the *.kcount files"
        "--file_names", "-f"
		help = "Text file of file name prefixes: '\$(file_names[i])-\$(kmer_size)mers.kcount'"
        "--sequence_dir", "-s"
		help = "Path to the individual FASTA sequence files"
        "--kmer_size", "-k"
		help = "Kmer size"
	"--jellyfish_location", "-j"
		help = "Location of the Jellyfish binary"
    end
    return parse_args(s)
end

#encode a kmer to an Int64, using 2 bits per character
@everywhere function encode(kmer)
	i = int64(0)
	bits = 0
	for c in kmer
       	if c == 'A'
			bits = 0
		elseif c == 'C'
			bits = 1
		elseif c == 'T'
			bits = 2
		elseif c == 'G'
			bits = 3
		else
            		error("Found non ACTG character in kmer")
		end
		i = (i << 2) + bits
	end
	return i
end

#Function that takes the jellyfish dump file and forms and writes the kcount file
@everywhere function kmer2kcount(kmers, output_dir, sequence_file, kmer_size)
	full_output_file = string(output_dir, "/", "$(sequence_file)-$(kmer_size)mers.kcount");
	all_lines_split = map(split,kmers); 
	kmers_to_write=zeros(Int64,(2,length(kmers)));
	kmer_words=Array(ASCIIString,length(kmers));
	kmer_words = [all_lines_split[i][1] for i=1:length(kmers)];
	counts = [all_lines_split[i][2] for i=1:length(kmers)];
	kmers_to_write[1,:] = map(encode,kmer_words);
	kmers_to_write[2,:] = map(int, counts);
	total = sum(kmers_to_write[2,:]);
	kmers_to_write = sortcols(kmers_to_write,alg=QuickSort) #takes about a minute...QuickSort is the fastest
	if isfile(full_output_file)
		rm(full_output_file)
	end
	fid = HDF5.h5open(full_output_file,"w")
	write(fid,"kmers",kmers_to_write)
	write(fid,"total_count",total)
	write(fid,"kmer_size",kmer_size)
	close(fid)
end


#This will form the jellyfish file, then dump it, then return just the kmers and counts themselves
@everywhere function fasta2kmer(sequence_dir, sequence_file, kmer_size, jellyfish_location)
	full_sequence_path = string(sequence_dir, "/", sequence_file);
	#count the kmers, then dump them
	kmers=readlines(`$(jellyfish_location) count $(full_sequence_path) -m $(kmer_size) -t 1 -s 100M -C -o /dev/fd/1 `  |> `$(jellyfish_location) dump /dev/fd/0 -c -t`);
	return kmers
end


function main()
	#Read in the args
	@everywhere parsed_args = remotecall_fetch(1,()->parse_commandline())
	println("Parsed args:")
	for (arg,val) in parsed_args
	    println("  $arg  =>  $val")
	end
	if parsed_args["output_dir"]==nothing
	    error("Output directory is required, use --help to see usage")
	end
	if parsed_args["file_names"]==nothing
	    error("Text file of file name prefixes is required, use --help to see usage")
	end
	if parsed_args["sequence_dir"]==nothing
	    error("Must provide the sequences directory, use --help to see usage")
	end
	if parsed_args["kmer_size"]==nothing
	    error("kmer_size must be provided, use --help to see usage")
	end
	#Use the command line arguments to populate variables
	@everywhere output_dir = parsed_args["output_dir"]
	@everywhere kmer_size = int(parsed_args["kmer_size"])
	@everywhere fid = open(parsed_args["file_names"],"r")
	@everywhere file_names = split(readall(fid))
	@everywhere close(fid)
	@everywhere sequence_dir = parsed_args["sequence_dir"]
	@everywhere jellyfish_location = parsed_args["jellyfish_location"]

	#Now do everything in parallel....
	pmap(i->kmer2kcount(fasta2kmer(sequence_dir, file_names[i], kmer_size, jellyfish_location), output_dir, file_names[i], kmer_size), 1:length(file_names), err_retry=true, err_stop=false)

end

main()