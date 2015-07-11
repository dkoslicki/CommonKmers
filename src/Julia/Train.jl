# ==============================================================================
# Train.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Forms the training data required for the CommonKmers method.
# ==============================================================================
# Call with something like: julia -p 3 Train.jl -i FullFileNames.txt -o ../Output/ -b /nfs1/Koslicki_Lab/koslickd/Bcalm/bcalm/./bcalm -r /data/temp/ -j /nfs1/Koslicki_Lab/koslickd/jellyfish-2.2.0/bin/jellyfish -c ./count_in_file
using HDF5
using ArgParse


#This function will call the old pmap on just a few of the specified workers. Note that pmap has been updated to allow this as an option, but I'm stuck with the older version for now.
function pmap2(f, lst, workers_to_use)
    np = length(workers_to_use)  # determine the number of processes available
    n = length(lst)
    results = cell(n)
    i = 1
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx > n
                            break
                        end
                        results[idx] = remotecall_fetch(p, f, lst[idx])
                    end
                end
            end
        end
    end
    results
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
		"--input_files", "-i"
			help = "File containing the full paths to the training genomes."
		"--output_folder", "-o"
			help = "Common Kmers training data path. Example: /path/to/CommonKmersData/Bcalms/"
		"--bcalm_binary", "-b"
			help = "Location of the bcalm_binary. eg ~/bin/./bcalm"
		"--ram_disk_location", "-r"
			help = "path to a ram disk (or SSD). Bcalm uses lots of disk IO, so this should be a fast IO device."
		"--jellyfish_binary", "-j"
			help = "Location of the jellyfish binary"
		"--count_in_file_binary", "-c"
			help = "Location of the count_in_file_binary. eg ~/bin/./count_in_file"
        "--size_of_chunk", "-s"
            help = "Matrix is computed in chunks, this determines the chunk size"
            arg_type = Int
            default = 500
        "--jellyfish_threads", "-t"
        	arg_type = Int
        	help = "Number of instances of jellyfish to use. Too many will cause disk thrashing"
        	default = 20
    end
    return parse_args(s)
end

#Parse the args
@everywhere parsed_args = remotecall_fetch(1,()->parse_commandline())
@everywhere input_files = parsed_args["input_files"]
@everywhere output_folder = parsed_args["output_folder"]
@everywhere bcalm_binary = parsed_args["bcalm_binary"]
@everywhere ram_disk_location = parsed_args["ram_disk_location"]
@everywhere jellyfish_binary = parsed_args["jellyfish_binary"]
@everywhere count_in_file_binary = parsed_args["count_in_file_binary"]
@everywhere chunk_size = parsed_args["size_of_chunk"]
@everywhere jellyfish_threads = parsed_args["jellyfish_threads"]

#Read in the file names to be processed
@everywhere fid = open(input_files,"r")
@everywhere file_names = map(x->basename(x),split(readall(fid)))
close(fid)
@everywhere num_files = length(file_names)
@everywhere fid = open(input_files,"r")
@everywhere full_file_names = split(readall(fid))
close(fid)

#Write the basenames to the training directory
fid = open("$(output_folder)/FileNames.txt","w")
for i=1:length(file_names)
	write(fid,"$(file_names[i])\n")
end



#Form the kmer counts
if isdir("$(output_folder)/Counts")
	rm("$(output_folder)/Counts", recursive=true)
end
mkdir("$(output_folder)/Counts")
@everywhere counts_directory = "$(output_folder)/Counts";

#Run at max jellyfish_threads jellyfish processes (to stop disks from thrashing)
worker_ids = workers();
if length(worker_ids)>jellyfish_threads
	workers_to_use = worker_ids[1:jellyfish_threads];
else
	workers_to_use = worker_ids;
end
pmap2(x->run(`$(jellyfish_binary) count $(full_file_names[x]) -m 30 -t 1 -s 100M -C -o $(output_folder)/Counts/$(file_names[x])-30mers.jf`),[1:num_files],workers_to_use);
pmap2(x->run(`$(jellyfish_binary) count $(full_file_names[x]) -m 50 -t 1 -s 100M -C -o $(output_folder)/Counts/$(file_names[x])-50mers.jf`),[1:num_files],workers_to_use);

#Form the CommonKmer matrix
#This will need to be refined after moving to count_in_file. Specifically, we should go row by row, and parallelize over columns. Due to the @sync, given a slow row, this can severely slow down the whole process
for kmer_size=[30;50]
	CommonKmersMatrix = SharedArray(Int64, (num_files,num_files), init=0);
	for j = 1:chunk_size:num_files
		for i = 1:chunk_size:num_files
			#only do the subdiagonal
			if i>=j
				#tic()
				#print("On (row,column) chunk: ($(i),$(j)) of $(num_files)\n")
				#parallelize over the chunk
				@sync begin
				@parallel for ii = i:i+chunk_size-1
					if ii<=num_files
						ikmers = "$(output_folder)/Counts/$(file_names[ii])-$(kmer_size)mers.jf"
						for jj=j:j+chunk_size-1
							if jj<=num_files
								jkmers = "$(output_folder)/Counts/$(file_names[jj])-$(kmer_size)mers.jf"
								(CommonKmersMatrix[ii,jj],CommonKmersMatrix[jj,ii]) = int64(split(readall(`$(count_in_file_binary) $(ikmers) $(jkmers)`)));
							end
						end
					end
				end
				end
				#toc()
				gc()
			end
		end
	end
	#Convert the matrix into something that HDF5 can write
	CommonKmersMatrix2=convert(Array,CommonKmersMatrix)'
	if isfile("$(output_folder)/CommonKmerMatrix-$(kmer_size)mers.h5")
		rm("$(output_folder)/CommonKmerMatrix-$(kmer_size)mers.h5")
	end
	#Write the output file
	h5write("$(output_folder)/CommonKmerMatrix-$(kmer_size)mers.h5","/common_kmers",CommonKmersMatrix2)
end

#Form Bcalm function
@everywhere function formBcalm(input_file, output_folder, bcalm_binary, ram_disk_location,jellyfish_binary,counts_directory)
	#Make temporary directory
	if isdir("$(ram_disk_location)/$(input_file)")
		rm("$(ram_disk_location)/$(input_file)",recursive=true)
	end
	mkdir("$(ram_disk_location)/$(input_file)")
	#Make .dot file
	run(`$(jellyfish_binary) dump $(counts_directory)/$(input_file)-30mers.jf -c -t` |> `cut -f 1` |> `tr '[:upper:]' '[:lower:]'` |> `sed 's/$/;/g'` |> "$(ram_disk_location)/$(input_file)/$(input_file)-30mers.dot")
	#Run Bcalm
	working_dir = pwd();
	cd("$(ram_disk_location)/$(input_file)")
	temp=readall(`$(bcalm_binary) $(input_file)-30mers.dot $(input_file)-30mers.bcalm 5`);
	#Move the file and delete the temp directory
	cd(working_dir)
	run(`cat $(ram_disk_location)/$(input_file)/$(input_file)-30mers.bcalm` |> `sed 's/;//g'` |> `tr '[:lower:]' '[:upper:]'` |> `sed 's/^/>seq\n/g'` |> "$(output_folder)/$(input_file)-30mers.bcalm.fa")
	rm("$(ram_disk_location)/$(input_file)",recursive=true)
end


#Run bcalm on all the files
if isdir("$(output_folder)/Bcalms")
	rm("$(output_folder)/Bcalms",recursive=true)
end
mkdir("$(output_folder)/Bcalms")
pmap(x->formBcalm(file_names[x], "$(output_folder)/Bcalms", bcalm_binary, ram_disk_location, jellyfish_binary, counts_directory),[1:num_files]);

#Remove Kmer counts
rm("$(output_folder)/Counts", recursive=true)


