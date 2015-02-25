#This is the function to classify the grinder experiments
#This uses the bcalm compressed kmer counts


#This assumes that the jellyfish djf files have already been formed, dumped, AND SORTED!
#jellyfish count SRR172902Even.fasta -m 30 -t 10 -s 100M -C -o SRR172902Even.fasta-30mers.jf
#jellyfish dump SRR172902Even.fasta-30mers.jf -c -t -o SRR172902Even.fasta-30mers.djf
#LC_ALL=C sort --buffer-size=10G --temporary-directory=/data/temp SRR172902Even.fasta-30mers.djf -o SRR172902Even.fasta-30mers.djf
#jellyfish count SRR172902Even.fasta -m 50 -t 10 -s 100M -C -o SRR172902Even.fasta-50mers.jf
#jellyfish dump SRR172902Even.fasta-50mers.jf -c -t -o SRR172902Even.fasta-50mers.djf
#I will first need to hash it, then sort it.
#First hash the 50mers djf file, then sort it, then form the the kcounts file
#HashDJFFunction(output_dir, sample_name, 50)
#LC_ALL=C sort --buffer-size=10G --temporary-directory=/data/temp SRR172902Even.fasta-50mers.djf -o SRR172902Even.fasta-50mers.djf

#Removed the stuff about sorting in linux


#Might be able to use the eval() method of forming djf file, as I'm only doing this on a few files, the number of file handles won't matter as much

#Load the required files
include("lsqnonneg.jl")
include("ConvertToCAMIOutputLCAFunctionDefault.jl")
include("ConvertToCAMIOutputLCAFunctionSensitive.jl")
include("ConvertToCAMIOutputLCAFunctionSpecific.jl")
using HDF5
using ArgParse

#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
		"--data_dir", "-d"
			help = "directory containing all the training files"
		"--output_dir", "-o"
			help = "directory where the output files will be written"
		"--input_file_name", "-i"
			help = "input file name"
		"--kind", "-k"
			help = "Type of output file. Options are: sensitive, specific, and default."
			default = "default"
		"--jellyfish_binary", "-j"
			help = "Location of jellyfish binary"
			default = "/home/pi/koslickd/jellyfish-2.1.1/bin/./jellyfish"
    end
    return parse_args(s)
end

#Parse the args
@everywhere parsed_args = remotecall_fetch(1,()->parse_commandline())
@everywhere data_dir = parsed_args["data_dir"]
@everywhere output_dir = parsed_args["output_dir"]
@everywhere input_file_name = parsed_args["input_file_name"]
@everywhere kind = parsed_args["kind"]
@everywhere jellyfish_binary = parsed_args["jellyfish_binary"]

#Set the input/output files
@everywhere file_names_path = "$(data_dir)/UniqueSpeciesFileNamesPruned.txt";
@everywhere taxonomy_file = "$(data_dir)/UniqueSpeciesTaxonomyPruned.txt";
@everywhere A30_file = "$(data_dir)/RepoPhlAn-12-20-14-UniqueSpeciesPruned-CommonKmerMatrix-30mersC.h5";
@everywhere A50_file = "$(data_dir)/RepoPhlAn-12-20-14-UniqueSpeciesPruned-CommonKmerMatrix-50mersC.h5";
@everywhere x_file = "$(output_dir)/$(basename(input_file_name))_reconstruction.txt"
@everywhere thresholds=[.90,.80,.70,.60,.50,.40,.30,.20,.10];
@everywhere normalize = "y";
@everywhere classification_file = "$(output_dir)/$(basename(input_file_name))_CommonKmers_classification.txt"
@everywhere num_threads = length(workers());

#form the jf files
run(`$(jellyfish_binary) count $(input_file_name) -m 30 -t $(num_threads) -s 100M -C -o $(output_dir)/$(basename(input_file_name))-30mers.jf`);
run(`$(jellyfish_binary) count $(input_file_name) -m 50 -t $(num_threads) -s 100M -C -o $(output_dir)/$(basename(input_file_name))-50mers.jf`);

#Form the Y functions
@everywhere fid = open(file_names_path,"r");
@everywhere file_names = split(readall(fid));
close(fid);
@everywhere num_files = length(file_names);
#do it once to read the jf and bcalms into memory
run(`$(data_dir)./query_per_sequence $(output_dir)/$(basename(input_file_name))-30mers.jf $(data_dir)/Bcalms/$(file_names[1])-30mers.bcalm.fa`);
Y30 = pmap(x->int(readall(`$(data_dir)./query_per_sequence $(output_dir)/$(basename(input_file_name))-30mers.jf $(data_dir)/Bcalms/$(file_names[x])-30mers.bcalm.fa`)),[1:num_files]);
run(`$(data_dir)./query_per_sequence $(output_dir)/$(basename(input_file_name))-50mers.jf $(data_dir)/Bcalms/$(file_names[1])-30mers.bcalm.fa`);
Y50 = pmap(x->int(readall(`$(data_dir)./query_per_sequence $(output_dir)/$(basename(input_file_name))-50mers.jf $(data_dir)/Bcalms/$(file_names[x])-30mers.bcalm.fa`)),[1:num_files]);
y30 = Y30/sum(Y30);
y50 = Y50/sum(Y50);

#Make the hypothetical matrices (later, only do this for the basis elements)
#30mers
A = float(h5read(A30_file,"/common_kmers"))'; #This was before I was transposing everything
A_norm = A./diag(A)';

#Create hypothetical organisms
hypothetical_matrix = zeros(size(A_norm,1),length(thresholds)*size(A_norm,2));
for threshold_index = 1:length(thresholds)
    threshold = thresholds[threshold_index];
    for organism_index = 1:size(A_norm,2)
        hyp_organism_vect = A_norm[:,organism_index]; #temp vector
        hyp_organism_vect[hyp_organism_vect .> threshold] = threshold; #round down to threshold
        hypothetical_matrix[:, organism_index + (size(A_norm,2)*(threshold_index-1))] = hyp_organism_vect; #set the corresponding column
    end
end
# real plus hypothetical
A_with_hypothetical30 = hcat(A_norm, hypothetical_matrix);
#50mers
A = float(h5read(A50_file,"/common_kmers"))'; 
A_norm = A./diag(A)';
#Create hypothetical organisms
thresholds=thresholds.^1.5;
hypothetical_matrix = zeros(size(A_norm,1),length(thresholds)*size(A_norm,2));
for threshold_index = 1:length(thresholds)
    threshold = thresholds[threshold_index];
    for organism_index = 1:size(A_norm,2)
        hyp_organism_vect = A_norm[:,organism_index]; #temp vector
        hyp_organism_vect[hyp_organism_vect .> threshold] = threshold; #round down to threshold
        hypothetical_matrix[:, organism_index + (size(A_norm,2)*(threshold_index-1))] = hyp_organism_vect; #set the corresponding column
    end
end
# real plus hypothetical
A_with_hypothetical50 = hcat(A_norm, hypothetical_matrix);
#together
A_with_hypothetical = vcat(A_with_hypothetical30, A_with_hypothetical50);


#perform the classification, non-restriction and non-sparsity promoting.
#set BLAS threads
blas_set_num_threads(length(workers()))
#tic();
#x = lsqnonneg(A_with_hypothetical, y);
#timing = toc();

#Perform the classification, sparsity promoting
lambda=1000;
basis=find(y30.>.001);
y = float(vcat(y30[basis],y50[basis]));
#tic();
x=lsqnonneg([ones(1,size(A_with_hypothetical,2));lambda*A_with_hypothetical[vcat(basis,basis.+num_files),:]],[0;lambda*y]);
#timing = toc();

#Normalize the result
x=x/sum(x);

#Write x_file
fid = open(x_file,"w")
for i=1:length(x)
	write(fid, "$(x[i])\n")
end
close(fid)

#Convert to CAMI format, use 30mer matrix to do the LCA
A = A_with_hypothetical[1:num_files, 1:num_files];
#A = float(h5read(Common_30mers_file,"/common_kmers"));
#A_norm = A./diag(A)';
sample_ID = "SAMPLEID";
contestant_ID = "CONTESTANTID";


if kind == "default"
	ConvertToCAMIOutputLCAFunctionDefault(x_file, taxonomy_file, classification_file, thresholds, A_norm, sample_ID, contestant_ID)
elseif kind == "specific"
	ConvertToCAMIOutputLCAFunctionSpecific(x_file, taxonomy_file, classification_file, thresholds, A_norm, sample_ID, contestant_ID)
elseif kind == "sensitive"
	ConvertToCAMIOutputLCAFunctionSensitive(x_file, taxonomy_file, classification_file, thresholds, A_norm, sample_ID, contestant_ID)
else
	error("Must choose one of the following output kinds for -t: default, specific, sensitive")
end


#Clean up the files
rm(x_file)
rm("$(output_dir)/$(basename(input_file_name))-30mers.jf")
rm("$(output_dir)/$(basename(input_file_name))-50mers.jf")








