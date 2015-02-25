# ==============================================================================
# ConvertToCAMIOutputLCA.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Takes a raw reconstruction vector (flat text file of floats on the same basis
# as the given taxonomy file) and outputs the classification in the CAMI format.
# ==============================================================================

#julia ConvertToCAMIOutputLCA.jl -i /nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAn/Samples/Classifications/testx_LCA.txt -t /nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAn/Taxonomy/FirstUniqueSpeciesFileNamesUniqueTaxonomyTaxID.txt -l 0 -o /nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAn/Samples/Classifications/test_classification_LCA.txt

function ConvertToCAMIOutputLCAFunctionSpecific(input_file, taxonomy_file, output_file, hyp_threshes, common_kmer_matrix_normalized, sample_ID, contestant_ID)

#Parse the args

#Constants, make these options next
#hyp_threshes = [.9, .8, .7, .6, .5, .4, .3, .2, .1]
cutoff = .001

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
#common_kmer_file = "/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAn/Output/FirstUniqueSpeciesFileNamesUnique-CommonKmers-40C.h5";
#common_kmer_matrix = float(h5read(common_kmer_file,"/common_kmers"))';
#common_kmer_matrix_normalized = common_kmer_matrix./diag(common_kmer_matrix)';

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
#		column = common_kmer_matrix_normalized[corresponding_real_organism_index,:];
		column = common_kmer_matrix_normalized[:,corresponding_real_organism_index];
		hyp_bin = int(floor((support_index-1)/num_organisms)); #This is the bin it belongs to #Make sure to add a -1 so that we use the right bin when the num_organisms evenly divides support_index
		hyp_thresh = hyp_threshes[hyp_bin]; #This is the corresponding threshold
	
		#LCA here
		temp = column;
#		temp[corresponding_real_organism_index] = 0; #set the organism entry to 0. Don't do this.
		distances = abs(temp.-hyp_thresh); #distances between column and thresh
		pos = indmin(distances); #This is the position of the nearest real organism to the threshold

		# If this distance is less than a bin's width, then do the LCA bit, otherwise go to a fixed taxonomy bit
		if distances[pos]<=0.2 #Later, make this adaptable to given input bins
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

end