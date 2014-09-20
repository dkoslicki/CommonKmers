#!/usr/local/julia/usr/bin/julia

###########################
#To Do
#1. Make sure that the abundances at each taxonomic rank sum to <=1
#2. Get it fully in the CAMI format (with TaxID's and ranks and TaxPath) see https://github.com/CAMI-challenge/contest_information/blob/master/file_formats/Example_TaxProfiling_Outputfile.txt
############################

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
			help = "Output evolutionary related level, the lth entry of [1, .95, .8, .7, .6, .5, .4, .3, .2, .1]"
			default = 1
        "--output_taxonomic_rank", "-r"
			help = "Output taxonomic rank, either an integer: 1, or a range [1:2]"
			default = "[1:8]"
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

#Just select the ones in the output level of interest
indicies_of_interest = range((output_level-1) * num_organisms, output_level * num_organisms)+1

#First, select the portion of the taxonomy that has a nonzero entry in the reconstruction
cutoff = .00001
support = indicies_of_interest[find(input[indicies_of_interest] .> cutoff)]
nonzero_taxonomy = taxonomy[(support % num_organisms).+1]

#open the output file
output_file_handle = open(output_file,"w")

#Write the header
write(output_file_handle,"# CAMI Submission for Taxonomic Profiling\n")
write(output_file_handle, "@Task:Profiling\n")
write(output_file_handle,"@Version:1.0\n")
write(output_file_handle,"@ContestantID:CONTESTANTID\n")
write(output_file_handle,"@SampleID:$(sample_ID)\n")
write(output_file_handle,"@Ranks: superkingdom|phylum|class|order|family|genus|species|strain\n")

#Now for each of the ranks, get the unique names, then loop over the non-zero taxonomy, increasing the value of the unique taxa name, then output these to the file
if typeof(output_taxonomic_rank) == Int64
	taxa_rank_list = [output_taxonomic_rank]
elseif typeof(output_taxonomic_rank) == Array{Int64,1}
	taxa_rank_list = output_taxonomic_rank
else
	error("Input taxonomic rank should be an integer, or else dont include the option to output all ranks")
end

for taxa_rank = taxa_rank_list	
	taxa_names = cell(0);
	for taxonomy_string = nonzero_taxonomy
		nonzero_taxonomy_split = split(taxonomy_string,"|")
		if length(nonzero_taxonomy_split) > taxa_rank
			taxa_name = join(nonzero_taxonomy_split[1:taxa_rank],"|")
			append!(taxa_names, {taxa_name})
		end
	end
	unique_taxa_names = sort(unique(taxa_names));
	#Now loop through each of the non_zero taxonomies, see if the taxa name matches, and then add this to the abundances
	taxa_abundances = Dict();
	for unique_taxa_name = unique_taxa_names
			taxa_abundances[unique_taxa_name] = 0
	end
	nonzero_taxonomy_counter = 1
	for taxonomy_string = nonzero_taxonomy
		nonzero_taxonomy_split = split(taxonomy_string,"|")
		if length(nonzero_taxonomy_split) > taxa_rank
			taxa_name = join(nonzero_taxonomy_split[1:taxa_rank],"|")
			taxa_abundances[taxa_name] = taxa_abundances[taxa_name] + input[support[nonzero_taxonomy_counter]]
		end
		nonzero_taxonomy_counter = nonzero_taxonomy_counter + 1
	end
	for unique_taxa_name = unique_taxa_names
		write(output_file_handle, "$(unique_taxa_name)")
		write(output_file_handle, "\t")
		write(output_file_handle, "$(taxa_abundances[unique_taxa_name])")
		write(output_file_handle, "\n")
	end
end

#Close the output file
close(output_file_handle)





