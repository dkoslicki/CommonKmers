# ==============================================================================
# ConvertToCAMIOutputLCA.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Takes a raw reconstruction vector (flat text file of floats on the same basis
# as the given taxonomy file) and outputs the classification in the CAMI format.
# ==============================================================================

#julia ConvertToCAMIOutputLCA.jl -i /nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAn/Samples/Classifications/testx_LCA.txt -t /nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAn/Taxonomy/FirstUniqueSpeciesFileNamesUniqueTaxonomyTaxID.txt -l 0 -o /nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAn/Samples/Classifications/test_classification_LCA.txt

##################################################
#Problems
#1. Things like 
#	10239	superkingdom	10239	Viruses	0.00014223541777389574
#	10239	phylum	10239|	Viruses|	3.950200757516657e-5
#Note the phylum and superkingdom have the same taxID
#2. Things like
#	1379698	strain	2||2||2||1379697|1379698	Bacteria||Bacteria||Bacteria||candidate_division_ZIXI|candidate_division_ZIXI_bacterium_RBG_1	3.589021327883041e-5


using HDF5
using ArgParse

#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input_file", "-i"
			help = "Input raw reconstruction text file"
        "--taxonomy_file", "-t"
			help = "taxonomy file, the ith line is the taxonomy of the ith training organism"
        "--output_level", "-l"
			help = "Output evolutionary related level, the lth entry of [1, .95, .8, .7, .6, .5, .4, .3, .2, .1]. Level of 0 means to sum over all of them."
			default = 0
        "--output_taxonomic_rank", "-r"
			help = "Output taxonomic rank, either an integer: 1, or a range [1:2]"
			default = "[1:7]"
		"--output_file", "-o"
			help = "Output text file"
		"--sample_ID", "-I"
			help = "Sample ID"
			default = "SAMPLEID"
		"--contestant_ID", "-C"
			help = "Contestant ID"
			default = "CONTESTANTID"
    end
    return parse_args(s)
end

#Parse the args
parsed_args = parse_commandline()
input_file = parsed_args["input_file"]
taxonomy_file = parsed_args["taxonomy_file"]
output_level = int(parsed_args["output_level"])
output_taxonomic_rank = eval(parse(parsed_args["output_taxonomic_rank"]))
output_file = parsed_args["output_file"]
sample_ID = parsed_args["sample_ID"]
contestant_ID = parsed_args["contestant_ID"]

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

##New

hyp_threshes = [.9, .8, .7, .6, .5, .4, .3, .2, .1]
cutoff = .00001
support = find(input .> cutoff)
#Create a massive taxonomy database
#Do the non-hypothetical organisms first, populating the species and strains
#Then do the the hypothetical organisms, doing the LCA, populating some of the higher taxonomy levels
#Then starting at the bottom, increment (or create, as the case may be) the higher taxonomic levels
#Then write the output.
common_kmer_file = "/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAn/Output/FirstUniqueSpeciesFileNamesUnique-CommonKmers-40C.h5";
common_kmer_matrix = float(h5read(common_kmer_file,"/common_kmers"))';
common_kmer_matrix_normalized = common_kmer_matrix./diag(common_kmer_matrix)';

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
		column = common_kmer_matrix_normalized[corresponding_real_organism_index,:];
		hyp_bin = int(floor(support_index/num_organisms)); #This is the bin it belongs to
		hyp_thresh = hyp_threshes[hyp_bin]; #This is the corresponding threshold
	
		#LCA here
		temp = column;
		temp[corresponding_real_organism_index] = 0; #set the organism entry to 0
		distances = abs(temp.-hyp_thresh); #distances between column and thresh
		pos = indmin(distances); #This is the position of the nearest real organism to the threshold
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
write(output_file_handle, "@Task:Profiling\n")
write(output_file_handle,"@Version:1.0\n")
write(output_file_handle,"@ContestantID:CONTESTANTID\n")
write(output_file_handle,"@SampleID:$(sample_ID)\n")
write(output_file_handle,"@Ranks: superkingdom|phylum|class|order|family|genus|species|strain\n")
write(output_file_handle,"\n")
write(output_file_handle,"@@TAXID\tRANK\tTAXPATH\tTAXPATH_SN\tPERCENTAGE\n")

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
		write(output_file_handle, "$(output_taxonomy[taxa_name])")
		write(output_file_handle, "\n")
	end
end

#Close the output file
close(output_file_handle)

exit()