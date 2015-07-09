# ==============================================================================
# FormBcalms.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Forms the Bcalm files required for the CommonKmers method.
# ==============================================================================
# call with something like: julia -p 2 FormBcalms.jl -i FileNames.txt -c ../Data/ -o ../Output/Bcalms/ -b /nfs1/Koslicki_Lab/koslickd/Bcalm/bcalm/./bcalm -r /nfs1/Koslicki_Lab/koslickd/CAMI/Training/Temp/ -j /nfs1/Koslicki_Lab/koslickd/jellyfish-2.2.0/bin/jellyfish
using ArgParse

@everywhere function formBcalm(input_file, output_folder, bcalm_binary, ram_disk_location,jelly_fish_binary,counts_directory)
	#Make temporary directory
	if isdir("$(ram_disk_location)/$(input_file)")
		rm("$(ram_disk_location)/$(input_file)",recursive=true)
	end
	mkdir("$(ram_disk_location)/$(input_file)")
	#Make .dot file
	run(`$(jelly_fish_binary) dump $(counts_directory)/$(input_file)-30mers.jf -c -t` |> `cut -f 1` |> `tr '[:upper:]' '[:lower:]'` |> `sed 's/$/;/g'` |> "$(ram_disk_location)/$(input_file)/$(input_file)-30mers.dot")
	#Run Bcalm
	working_dir = pwd();
	cd("$(ram_disk_location)/$(input_file)")
	temp=readall(`$(bcalm_binary) $(input_file)-30mers.dot $(input_file)-30mers.bcalm 5`);
	#Move the file and delete the temp directory
	cd(working_dir)
	run(`cat $(ram_disk_location)/$(input_file)/$(input_file)-30mers.bcalm` |> `sed 's/;//g'` |> `tr '[:lower:]' '[:upper:]'` |> `sed 's/^/>seq\n/g'` |> "$(output_folder)/$(input_file)-30mers.bcalm.fa")
	rm("$(ram_disk_location)/$(input_file)",recursive=true)
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
		"--input_files", "-i"
			help = "File containing base names of training genomes."
		"--counts_directory", "-c"
			help = "Directory containing 30mer counts. Files MUST be of the for file_name-30mers.jf"
		"--output_folder", "-o"
			help = "Common Kmers training data path. Example: /path/to/CommonKmersData/Bcalms/"
		"--bcalm_binary", "-b"
			help = "Location of the bcalm_binary. eg ~/bin/./bcalm"
		"--ram_disk_location", "-r"
			help = "path to a ram disk (or SSD)."
		"--jelly_fish_binary", "-j"
			help = "Location of the jellyfish binary"
    end
    return parse_args(s)
end

#Parse the args
@everywhere parsed_args = remotecall_fetch(1,()->parse_commandline())
@everywhere input_files = parsed_args["input_files"]
@everywhere output_folder = parsed_args["output_folder"]
@everywhere bcalm_binary = parsed_args["bcalm_binary"]
@everywhere ram_disk_location = parsed_args["ram_disk_location"]
@everywhere jelly_fish_binary = parsed_args["jelly_fish_binary"]
@everywhere counts_directory = parsed_args["counts_directory"]

#Check for errors
if parsed_args["input_files"]==nothing
    error("A file containing the full paths to the training genomes is required, use --help to see usage")
end
if parsed_args["output_folder"]==nothing
    error("An output folder is required, use --help to see usage")
end
if parsed_args["bcalm_binary"]==nothing
    error("bcalm_binary location is required, use --help to see usage")
end
if parsed_args["ram_disk_location"]==nothing
    error("ram_disk_location path is required, use --help to see usage")
end
if parsed_args["jelly_fish_binary"]==nothing
    error("jelly_fish_binary location is required, use --help to see usage")
end
if parsed_args["counts_directory"]==nothing
    error("counts_directory path is required, use --help to see usage")
end

#Read in the file names to be processed
@everywhere fid = open(input_files,"r")
@everywhere file_names = split(readall(fid))
close(fid)
@everywhere num_files = length(file_names)

#Run bcalm on all the files
pmap(x->formBcalm(file_names[x], output_folder, bcalm_binary, ram_disk_location,jelly_fish_binary,counts_directory),[1:num_files]);









