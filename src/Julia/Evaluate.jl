# ==============================================================================
# Evaluate.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Creates a file of error that outputs various error metrics at various taxonomic
# ranks for a test file and a "ground truth" file.
# ==============================================================================

using ArgParse

#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input_file", "-i"
			help = "Input CAMI reconstruction file"
        "--ground_truth", "-g"
			help = "Input CAMI file for ground truth"
		"--output_file", "-o"
			help = "Output text file"
		"--normalize", "-n"
			help = "y/n to normalize the inputs"
			default = "n"
		"--epsilon", "-e"
			help = "Cutoff for determining if organism is actually present in sample. Used for all metrics. Default is 0"
			default = "0"
		"--num_classes", "-c"
			help = "Number of possible classification classes, this is needed to calculate the True Negative rate. If you can classify to X different organisms, then us -c X"
    end
    return parse_args(s)
end



#Parse the args
parsed_args = parse_commandline()
input_file = parsed_args["input_file"]
ground_truth = parsed_args["ground_truth"]
output_file = parsed_args["output_file"]
normalize = parsed_args["normalize"]
epsilon = float(parsed_args["epsilon"])
num_classes = int(parsed_args["num_classes"])

#Read input
fid = open(input_file,"r")
input_split = map(x->split(x),readlines(fid));
close(fid)
input_len = length(input_split)
#Get rid of the header
iter = 1;
while iter <= input_len
	if length(input_split[iter]) >= 1
		if input_split[iter][1] == "@@TAXID"
			break
		end
	end
	iter = iter + 1
end
if iter >= input_len
	error("Missing @@TAXID from file $(input_file)\n")
end

iter = iter + 1;
#TaxID's
input_taxIDs = int(map(x->input_split[x][1],iter:input_len))
#Ranks
input_ranks = map(x->input_split[x][2],iter:input_len)
#Frequencies
input_freqs = float(map(x->input_split[x][end],iter:input_len))

#Read the ground truth
fid = open(ground_truth,"r")
gt_split = map(x->split(x),readlines(fid));
close(fid)
gt_input_len = length(gt_split)
#Get rid of the header
iter = 1;
while iter <= gt_input_len
	if length(gt_split[iter]) >= 1
		if gt_split[iter][1] == "@@TAXID"
			break
		end
	end
	iter = iter + 1
end
if iter >= input_len
	error("Missing @@TAXID from file: $(ground_truth)\n")
end

iter = iter + 1;
#TaxID's
gt_taxIDs = int(map(x->gt_split[x][1],iter:gt_input_len))
#Ranks
gt_ranks = map(x->gt_split[x][2],iter:gt_input_len)
#Frequencies
gt_freqs = float(map(x->gt_split[x][end],iter:gt_input_len))

#Use cutoff
input_freqs[input_freqs .< epsilon] = 0.
gt_freqs[gt_freqs .< epsilon] = 0.

#Get the unique ranks
unique_rank_names = unique(vcat(gt_ranks,input_ranks))

#open output file
fid = open(output_file,"w")
write(fid, "@input\t $(basename(input_file))\n")
write(fid, "@ground_truth\t $(basename(ground_truth))\n")
write(fid, "@Taxanomic ranks\t $(join(unique_rank_names,"|"))\n")

Precs = zeros(length(unique_rank_names))
Senss = zeros(length(unique_rank_names))
Specs = zeros(length(unique_rank_names))
Accs = zeros(length(unique_rank_names))
FPRs = zeros(length(unique_rank_names))
FDRs = zeros(length(unique_rank_names))
L1norms = zeros(length(unique_rank_names))
L2norms = zeros(length(unique_rank_names))

iter = 1
for rank in unique_rank_names
	#Select Tax ID's of interest from input
	input_rank_taxIDs = Int64[];
	input_rank_freqs = Float64[];
	for i=1:length(input_ranks)
		if input_ranks[i] == rank
			push!(input_rank_taxIDs, input_taxIDs[i])
			push!(input_rank_freqs, input_freqs[i])
		end
	end
	
	#Select Tax ID's of interest from ground truth
	gt_rank_taxIDs = Int64[];
	gt_rank_freqs = Float64[];
	for i=1:length(gt_ranks)
		if gt_ranks[i] == rank
			push!(gt_rank_taxIDs, gt_taxIDs[i])
			push!(gt_rank_freqs, gt_freqs[i])
		end
	end
	
	#Normalize if this was asked for
	if normalize == "y"
		input_rank_freqs = input_rank_freqs/sum(input_rank_freqs)
		gt_rank_freqs = gt_rank_freqs/sum(gt_rank_freqs)
	elseif normalize == "n"
		input_rank_freqs = input_rank_freqs
		gt_rank_freqs = gt_rank_freqs
	else
		error("--normalize option should be one of 'y' or 'n'")
	end
	
	#Calculate the errors
	set_gt_rank_taxIDs = Set(gt_rank_taxIDs)
	set_input_rank_taxIDs = Set(input_rank_taxIDs)
	#True positive
	TP = length(intersect(set_gt_rank_taxIDs, set_input_rank_taxIDs))
	FP = length(setdiff(set_input_rank_taxIDs, set_gt_rank_taxIDs))
	TN = num_classes - length(union(set_gt_rank_taxIDs, set_input_rank_taxIDs))
	FN = length(setdiff(set_gt_rank_taxIDs, set_input_rank_taxIDs))
	Prec = TP/(TP+FP)
	Sens = TP/(TP+FN)
	Spec = TN/(TN+FP)
	Acc = (TP+TN)/(TP+TN+FP+FN)
	FPR = 1 - Spec
	FDR = 1 - Prec
	
	#Populate output vector
	Precs[iter] = Prec
	Senss[iter] = Sens
	Specs[iter] = Spec
	Accs[iter] = Acc
	FPRs[iter] = FPR
	FDRs[iter] = FDR
	
	#L1 and L2 error
	all_taxIDs = unique(vcat(gt_rank_taxIDs,input_rank_taxIDs))
	L1sum = 0.
	L2sum = 0.
	#Currently this is very lazily (=inefficiently) done. What I should is pre-split the three cases, then do the sum. Whatever, this is fast enough
	for taxID = all_taxIDs
			#note that taxID's might repeat in gt_rank_taxIDs due to missing organism, so sum up all of these before taking the difference
			gt_sum = sum(gt_rank_freqs[findin(gt_rank_taxIDs,taxID)])
			input_sum = sum(input_rank_freqs[findin(input_rank_taxIDs,taxID)])
			L1sum = L1sum + abs(gt_sum - input_sum)
			L2sum = L2sum + abs(gt_sum - input_sum).^2
	end
	L1norm = L1sum
	L2norm = sqrt(L1sum)
	
	#Populate output vector
	L1norms[iter] = L1norm
	L2norms[iter] = L2norm
	
	iter = iter + 1;
end

#write to output file
write(fid, "Precision\t $(join(Precs,"|"))\n")
write(fid, "Sensitivity\t $(join(Senss,"|"))\n")
write(fid, "Specificity\t $(join(Specs,"|"))\n")
write(fid, "Accuracy\t $(join(Accs,"|"))\n")
write(fid, "False Positive Rate\t $(join(FPRs,"|"))\n")
write(fid, "False Discovery Rate\t $(join(FDRs,"|"))\n")
write(fid, "L1 norm\t $(join(L1norms,"|"))\n")
write(fid, "L2 norm\t $(join(L2norms,"|"))\n")

close(fid)
