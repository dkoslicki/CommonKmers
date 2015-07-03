using HDF5
using ArgParse


function lsqnonneg(C::Matrix, d::Vector, tol::Real=-1, itmax_factor::Real=3)
	#set the tolerance
	(m,n) = size(C);
	tol = (tol == -1) ? tol = 10*eps()*norm(C,1)*(maximum(size(C))+1) : tol
	itmax = itmax_factor*n;
	
	# Initialize vector of n zeros and Infs (to be used later)
	wz = zeros(n,1);
	
	# Initialize set of non-active columns to null
	P = falses(n,1);
	
	# Initialize set of active columns to all and the initial point to zeros
	Z = trues(n,1);
	x = zeros(n,1);
    Ctrans=C';
    dtemp=d[:];
    BLAS.gemv!('N',-1.,C,x[:],1.,dtemp); #resid = d - C*x;
    w=BLAS.gemv('N',1.,Ctrans, dtemp); #w = Ctrans*resid;
	
	# Set up iteration criterion
	outeriter = 0;
	iter = 0;
	exitflag = 1;
	
	# Outer loop to put variables into set to hold positive coefficients
	while any(Z) & any(w[Z[:]] .> tol)
		#print("On iteration $(outeriter)\n")

		outeriter = outeriter + 1; 
		
		# Reset intermediate solution z
		z = zeros(n,1);
		
		# Create wz, a Lagrange multiplier vector of variables in the zero set.
		# wz must have the same size as w to preserve the correct indices, so
		# set multipliers to -Inf for variables outside of the zero set.
		wz[P] = -Inf;
		wz[Z] = w[Z[:]];

		# Find variable with largest Lagrange multiplier
		t = indmax(wz);
		
		# Move variable t from zero set to positive set
		P[t] = true;
		Z[t] = false;
		
		# Compute intermediate solution using only variables in positive set
		z[P] = C[:,find(P)]\d;
		
		#inner loop to remove elements from the positive set which no longer belong
		while any(z[P] .<= 0)
			#print("entering inner loop\n")
			iter = iter + 1;
			if iter > itmax
				print("lsqnonneg:IterationCountExceeded");
		              exitflag = 0;
				iterations = outeriter;
		              resnorm = sum(resid.*resid);
				x = z;
				lambda = w;
				return x
			end
			
			# Find indices where intermediate solution z is approximately negative
			Q = [(z .<= 0) & P];

			# Choose new x subject to keeping new x nonnegative
			alpha = minimum(x[Q]./(x[Q] - z[Q]));
			x = x + alpha*(z - x);

			# Reset Z and P given intermediate values of x
			Z = [((abs(x) .< tol) & P) | Z ];
			P = ~Z;
			z = zeros(n,1);           # Reset z
			z[P] = C[:,find(P)]\d;      # Re-solve for z
		end

		x = z;

        dtemp=d[:];
        BLAS.gemv!('N',-1.,C,x[:],1.,dtemp); #resid = d - C*x;
        w=BLAS.gemv('N',1.,Ctrans, dtemp); #w = Ctrans*resid;

	end
	return x
end


# ==============================================================================
# ConvertToCAMIOutputLCA.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Takes a raw reconstruction vector (flat text file of floats on the same basis
# as the given taxonomy file) and outputs the classification in the CAMI format.
# ==============================================================================

function ConvertToCAMIOutputLCA(kind, input_file, taxonomy_file, output_file, hyp_threshes, common_kmer_matrix_normalized, sample_ID, contestant_ID)

#Parse the args

#Constants, make these options next
#hyp_threshes = [.9, .8, .7, .6, .5, .4, .3, .2, .1]
if kind=="sensitive"
	cutoff = .0001
else
	cutoff = .001
end

#Read in the input file
fid = open(input_file,"r")
input = readlines(fid)
close(fid)
input = map(x->float(strip(x)), input)

#Next, read in the taxonomy file
fid = open(taxonomy_file,"r")
taxonomy = readlines(fid)
close(fid)
taxonomy = map(x->strip(split(x)[3]), taxonomy)
num_organisms = length(taxonomy)

#New
support = find(input .> cutoff)
#Create a massive taxonomy database
#Do the non-hypothetical organisms first, populating the species and strains
#Then do the the hypothetical organisms, doing the LCA, populating some of the higher taxonomy levels
#Then starting at the bottom, increment (or create, as the case may be) the higher taxonomic levels
#Then write the output.

output_taxonomy = Dict()
for support_index = support
	if support_index <= length(taxonomy) #If it's not a hypothetical organism
		if haskey(output_taxonomy, taxonomy[support_index]) #If it's in there, add to it.
			output_taxonomy[taxonomy[support_index]] = output_taxonomy[taxonomy[support_index]] + input[support_index]
		else
			output_taxonomy[taxonomy[support_index]] = input[support_index]
		end
	else #it's a hypothetical organism, so do the LCA here
		corresponding_real_organism_index = mod(support_index,num_organisms);
		#fix zero guy
		if corresponding_real_organism_index == 0
			corresponding_real_organism_index = num_organisms
		end
		column = common_kmer_matrix_normalized[:,corresponding_real_organism_index];
		hyp_bin = int(floor((support_index-1)/num_organisms)); #This is the bin it belongs to #Make sure to add a -1 so that we use the right bin when the num_organisms evenly divides support_index
		hyp_thresh = hyp_threshes[hyp_bin]; #This is the corresponding threshold
	
		#LCA here
		temp = column;
		distances = abs(temp.-hyp_thresh); #distances between column and thresh
		pos = indmin(distances); #This is the position of the nearest real organism to the threshold

		# If this distance is less than a bin's width, then do the LCA bit, otherwise go to a fixed taxonomy bit
		if kind=="default" || distances[pos]<=0.2 #Later, make this adaptable to given input bins
			#Now find the LCA between this organism and the hypothetical one
			hyp_taxonomy_split = split(taxonomy[corresponding_real_organism_index],"|"); #Hyp taxonomy based on corresponding organism
			candidate_LCA_taxnonmy_split = split(taxonomy[pos],"|"); #closest organism taxnonmy
			LCA = 1; #In case it only agrees to the kingdom level, and we need to declare LCA as a global variable
			for LCA_index = minimum([length(candidate_LCA_taxnonmy_split), length(hyp_taxonomy_split)]):-1:1
				if hyp_taxonomy_split[LCA_index] == candidate_LCA_taxnonmy_split[LCA_index]
					LCA = LCA_index
					break
				end
			end
		else #Fixed LCA based on bin
			hyp_taxonomy_split = split(taxonomy[corresponding_real_organism_index],"|");
			candidate_LCA_taxnonmy_split = hyp_taxonomy_split;
			LCA = 1;
			if kind=="sensitive"
				if hyp_thresh == .9
					LCA = length(hyp_taxonomy_split);
				elseif hyp_thresh == .8
					LCA = length(hyp_taxonomy_split);
				elseif hyp_thresh == .7
					LCA = maximum([1, length(hyp_taxonomy_split)-1]);
				elseif hyp_thresh == .6
					LCA = maximum([1, length(hyp_taxonomy_split)-1]);
				elseif hyp_thresh == .5
					LCA = maximum([1, length(hyp_taxonomy_split)-2]);
				elseif hyp_thresh == .4
					LCA = maximum([1, length(hyp_taxonomy_split)-2]);
				elseif hyp_thresh == .3
					LCA = maximum([1, length(hyp_taxonomy_split)-2]);
				elseif hyp_thresh == .2
					LCA = maximum([1, length(hyp_taxonomy_split)-3]);
				elseif hyp_thresh == .1
					LCA = maximum([1, length(hyp_taxonomy_split)-4]);
				end
			elseif kind=="specific"
				if hyp_thresh == .9
					LCA = length(hyp_taxonomy_split);
				elseif hyp_thresh == .8
					LCA = length(hyp_taxonomy_split);
				elseif hyp_thresh == .7
					LCA = maximum([1, length(hyp_taxonomy_split)-1]);
				elseif hyp_thresh == .6
					LCA = maximum([1, length(hyp_taxonomy_split)-2]);
				elseif hyp_thresh == .5
					LCA = maximum([1, length(hyp_taxonomy_split)-2]);
				elseif hyp_thresh == .4
					LCA = maximum([1, length(hyp_taxonomy_split)-2]);
				elseif hyp_thresh == .3
					LCA = maximum([1, length(hyp_taxonomy_split)-2]);
				elseif hyp_thresh == .2
					LCA = maximum([1, length(hyp_taxonomy_split)-3]);
				elseif hyp_thresh == .1
					LCA = maximum([1, length(hyp_taxonomy_split)-4]);
				end
			else
				error("invalid kind (-k) option. Options are: default, sensitive, specific")
			end
		end

		#Now update the taxonomy dictionary
		LCA_taxonomy = join(candidate_LCA_taxnonmy_split[1:LCA],"|")
		if haskey(output_taxonomy,LCA_taxonomy) #If it's in there, add to it.
			output_taxonomy[LCA_taxonomy] = output_taxonomy[LCA_taxonomy] + input[support_index]
		else
			output_taxonomy[LCA_taxonomy] = input[support_index]
		end
	end
end

#Now sum up from coarser taxonomic ranks to higher ones (summing finer to coarser for each taxa level)
temp_dict = copy(output_taxonomy);
for taxonomic_level = 2:8
	for key in keys(temp_dict)# or use: key = [key for key in keys(output_taxonomy)]. Can't use key in keys(output_taxonomy) because those keys change
		split_key = split(key,"|");
		if length(split_key) == taxonomic_level
			#loop through the higher taxonomic levels
			for higher_level = (taxonomic_level-1):-1:1
				higher_taxonomy = join(split_key[1:higher_level],"|")
				#update the dictionary value
				if haskey(output_taxonomy,higher_taxonomy) #If it's in there, add to it.
					output_taxonomy[higher_taxonomy] = output_taxonomy[higher_taxonomy] + output_taxonomy[key] #add the value at the base taxonomy
				else #If note, make it equal to the base taxonomy value
					output_taxonomy[higher_taxonomy] = output_taxonomy[key]
				end
			end
		end
	end
end


#open the output file
output_file_handle = open(output_file,"w")

#Write the header
write(output_file_handle,"# CAMI Submission for Taxonomic Profiling\n")
#write(output_file_handle, "@Task:Profiling\n")
write(output_file_handle,"@Version:0.9.1\n")
#write(output_file_handle,"@ContestantID:CONTESTANTID\n")
write(output_file_handle,"@SampleID:$(sample_ID)\n")
write(output_file_handle,"@Ranks: superkingdom|phylum|class|order|family|genus|species|strain\n")
write(output_file_handle,"\n")
write(output_file_handle,"@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")

taxa_names = [key for key in keys(output_taxonomy)];
taxa_names = sort(taxa_names, by=x->length(split(x,"|")));

for taxa_name = taxa_names
	taxID = split(split(taxa_name,"|")[end],"_")[3];
	
	rankAbvr = split(split(taxa_name,"|")[end],"_")[1];
	if rankAbvr == "k"
		rank = "superkingdom"
	elseif rankAbvr == "p"
		rank = "phylum"
	elseif rankAbvr == "c"
		rank = "class"
	elseif rankAbvr == "o"
		rank = "order"
	elseif rankAbvr == "f"
		rank = "family"
	elseif rankAbvr == "g"
		rank = "genus"
	elseif rankAbvr == "s"
		rank = "species"
	elseif rankAbvr == "t"
		rank = "strain"
	else
		rank = "unknown"
	end

		
	taxPath = map(x->split(x,"_")[3],split(taxa_name,"|")); #Tax ID's
	taxPathSN = map(x->join(split(x,"_")[4:end],"_"), split(taxa_name,"|")); #Taxa names
	#If a Tax ID is repeated at a lower taxonomic rank, this means that that rank is missing, so let's just delete it.
	for i=length(taxPath):-1:2
		if i>=2
			if taxPath[i] == taxPath[i-1]
				taxPath[i] = ""
				taxPathSN[i] = ""
			end
		end
	end
	#If it ends with a blank, we don't want to be including it
	if taxPath[length(taxPath)] != ""
		write(output_file_handle, "$(taxID)")
		write(output_file_handle, "\t")
		write(output_file_handle, "$(rank)")
		write(output_file_handle, "\t")
		#Join back up the paths
		taxPath = join(taxPath,"|")
		taxPathSN = join(taxPathSN,"|")
	
		write(output_file_handle, "$(taxPath)")
		write(output_file_handle, "\t")
	
		write(output_file_handle, "$(taxPathSN)")
		write(output_file_handle, "\t")
		write(output_file_handle, "$(100*output_taxonomy[taxa_name])") #Turn into a frequency (i.e. 90 not .9)
		write(output_file_handle, "\n")
	end
end

#Close the output file
close(output_file_handle)

end




#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
		"--data_dir", "-d"
			help = "directory containing all the training files"
		"--output_file", "-o"
			help = "output file"
		"--input_file_name", "-i"
			help = "input file name"
		"--kind", "-k"
			help = "Type of output file. Options are: sensitive, specific, and default."
			default = "default"
		"--jellyfish_binary", "-j"
			help = "Location of jellyfish binary. eg ~/bin/./jellyfish"
			default = "jellyfish"
		"--query_per_sequence_binary", "-q"
			help = "Location of the query_per_sequence binary. eg ~/bin/./query_per_sequence"
			default = "query_per_sequence"
		"--contestant_ID", "-c"
			help = "Contestant ID"
			default = "CONTESTANTID"
		"--quality", "-Q"
			help = "minimum per-base quality score required to include kmer"
			default = "C"
    end
    return parse_args(s)
end

#Parse the args
@everywhere parsed_args = remotecall_fetch(1,()->parse_commandline())
@everywhere data_dir = parsed_args["data_dir"]
@everywhere output_file = parsed_args["output_file"]
@everywhere input_file_name = parsed_args["input_file_name"]
@everywhere kind = parsed_args["kind"]
@everywhere jellyfish_binary = parsed_args["jellyfish_binary"]
@everywhere query_per_sequence_binary = parsed_args["query_per_sequence_binary"]
@everywhere contestant_ID = parsed_args["contestant_ID"];
@everywhere sample_ID = input_file_name
@everywhere quality = parsed_args["quality"]

#Set the input/output files
@everywhere file_names_path = "$(data_dir)/UniqueSpeciesFileNamesPruned.txt";
@everywhere taxonomy_file = "$(data_dir)/UniqueSpeciesTaxonomyPruned.txt";
@everywhere A30_file = "$(data_dir)/RepoPhlAn-12-20-14-UniqueSpeciesPruned-CommonKmerMatrix-30mersC.h5";
@everywhere A50_file = "$(data_dir)/RepoPhlAn-12-20-14-UniqueSpeciesPruned-CommonKmerMatrix-50mersC.h5";
@everywhere x_file = "$(basename(input_file_name))_reconstruction.txt"
@everywhere thresholds=[.90,.80,.70,.60,.50,.40,.30,.20,.10];
@everywhere normalize = "y";
#@everywhere classification_file = "$(output_dir)/$(basename(input_file_name))_CommonKmers_classification.txt"
@everywhere classification_file = output_file;
@everywhere num_threads = length(workers());

#form the jf files
run(`$(jellyfish_binary) count $(input_file_name) -m 30 -t $(num_threads) -s 100M -C -Q $(quality) -o $(basename(input_file_name))-30mers.jf`);
run(`$(jellyfish_binary) count $(input_file_name) -m 50 -t $(num_threads) -s 100M -C -Q $(quality) -o $(basename(input_file_name))-50mers.jf`);

#Form the Y functions
@everywhere fid = open(file_names_path,"r");
@everywhere file_names = split(readall(fid));
close(fid);
@everywhere num_files = length(file_names);
#do it once to read the jf and bcalms into memory
temp=readall(`$(query_per_sequence_binary) $(basename(input_file_name))-30mers.jf $(data_dir)/Bcalms50/$(file_names[1])-50mers.bcalm.fa`);
Y30 = pmap(x->int(readall(`$(query_per_sequence_binary) $(basename(input_file_name))-30mers.jf $(data_dir)/Bcalms50/$(file_names[x])-50mers.bcalm.fa`)),[1:num_files]);
#now for the 50mers
temp=readall(`$(query_per_sequence_binary) $(basename(input_file_name))-50mers.jf $(data_dir)/Bcalms50/$(file_names[1])-50mers.bcalm.fa`);
Y50 = pmap(x->int(readall(`$(query_per_sequence_binary) $(basename(input_file_name))-50mers.jf $(data_dir)/Bcalms50/$(file_names[x])-50mers.bcalm.fa`)),[1:num_files]);
y30 = Y30/float(split(readall(`$(jellyfish_binary) stats $(basename(input_file_name))-30mers.jf`))[6]); #divide by total number of kmers in sample
y50 = Y50/float(split(readall(`$(jellyfish_binary) stats $(basename(input_file_name))-50mers.jf`))[6]);

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

#Perform the classification, just lsqnonneg
y = float(vcat(y30,y50));
#tic();
x=lsqnonneg(A_with_hypothetical,y,.0005,3);
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
A = float(h5read(A30_file,"/common_kmers"))'; #This was before I was transposing everything
A_norm = A./diag(A)';


if kind=="default" || kind=="specific" || kind=="sensitive"
	ConvertToCAMIOutputLCA(kind, x_file, taxonomy_file, classification_file, thresholds, A_norm, sample_ID, contestant_ID)
else
	error("Must choose one of the following output kinds for -k: default, specific, sensitive")
end


#Clean up the files
rm(x_file)
rm("$(basename(input_file_name))-30mers.jf")
rm("$(basename(input_file_name))-50mers.jf")








