#This will plot the bubble chart on the tree
import random, math
from ete2 import Tree, faces, TreeStyle, COLOR_SCHEMES, TextFace, BarChartFace, CircleFace, AttrFace, Phyloxml, phyloxml
import sys, getopt

def main(argv):
	input_file=''
	title='Title'
	label_internal_nodes = False
	label_leaves = False
	out_file=''
	width=750
	out_file_xml=''
	try:
		opts, args = getopt.getopt(argv,"h:i:t:lno:w:x:",["Help=","InputFile=","Title=","LabelLeaves=", "LabelInternalNodes=","OutFile=","Width=","OutFileXML="])
	except getopt.GetoptError:
		print 'Unknown option, call using: ./PlotTree.py -i <InputCAMIFile> -t <Title> -l <LabelLeavesFlag> -n <LabelInternalNodesFlag> -o <OutFile.png> -x <Outfile.xml> -w <Width>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './PlotTree.py -i <InputCAMIFile> -t <Title> -l <LabelLeavesFlag> -n <LabelInternalNodesFlag> -o <OutFile> -x <OutFile.xml> -w <Width>'
			sys.exit(2)
		elif opt in ("-i", "--InputFile"):
			input_file = arg
		elif opt in ("-t", "--Title"):
			title = arg
		elif opt in ("-l", "--LabelLeaves"):
			label_leaves = True
		elif opt in ("-n","--LabelInternalNodes"):
			label_internal_nodes = True
		elif opt in ("-o", "--OutFile"):
			out_file = arg
		elif opt in ("-w", "--Width"):
			width = int(arg)
		elif opt in ("-x", "--OutFileXML"):
			out_file_xml = arg
			
	schema_names = COLOR_SCHEMES.keys()
	
	#Read the common kmer profile
	ckm_tax_paths = []
	ckm_name_to_perc = dict()
	fid = open(input_file,'r')
	file = fid.readlines()
	fid.close()
	
	#Put placeholders in for missing names like: "||" -> "|NA1|"
	file_noblank = list()
	i=0
	for line in file:
		while "||" in line:
			line = line.replace("||","|NONAME|",1)
			i = i+1
		file_noblank.append(line)
	
	#Get the names and weights
	for line in file_noblank:
		if line[0]!='#' and line[0]!='@' and line[0]!='\n': #Don't parse comments or blank lines
			temp = line.split()[3] #Get the names
			ckm_tax_paths.append(temp)
			ckm_name_to_perc[temp.split("|")[-1]] = line.split()[-1] #Get the weights
	
	#Create the tree
	t=Tree()
	names_to_nodes = dict()
	for i in range(0,len(ckm_tax_paths)):
		split_tax_path = ckm_tax_paths[i].split("|")
		if len(split_tax_path)==1: #If len==1, then it's a superkingdom
			names_to_nodes[split_tax_path[0]] = t.add_child(name=split_tax_path[0]) #connect directly to tree
		else:
			if split_tax_path[-2] in names_to_nodes: #If the parent is already in the tree, add to tree
				names_to_nodes[split_tax_path[-1]] = names_to_nodes[split_tax_path[-2]].add_child(name=split_tax_path[-1])
			else: #Otherwise iterate up until we have something that is in the tree
				j=2
				while split_tax_path[-j]=="NONAME":
					j = j + 1
				#This skips over the NONAMES
				names_to_nodes[split_tax_path[-1]] = names_to_nodes[split_tax_path[-j]].add_child(name=split_tax_path[-1])
	
	#Show the tree
	#print t.get_ascii(show_internal=True)
	
	#scheme = random.sample(schema_names, 1)[0] #'set2' is nice, 
	scheme = 'set2'

	def layout(node):
		if node.name in ckm_name_to_perc:
			ckm_perc = float(ckm_name_to_perc[node.name])
		else:
			ckm_perc = 0
		F = CircleFace(radius=3.14*math.sqrt(ckm_perc), color="RoyalBlue", style="sphere")
		F.border.width = None
		F.opacity = 0.6
		faces.add_face_to_node(F,node, 0, position="branch-right")
		if label_internal_nodes:
			faces.add_face_to_node(TextFace(node.name, fsize=7),node, 0, position="branch-top")
	
	ts = TreeStyle()
	ts.layout_fn = layout
	ts.mode = "r"
	ts.show_leaf_name = label_leaves
	ts.min_leaf_separation = 50
	ts.title.add_face(TextFace(title, fsize=20), column=0)
	
	#Export the tree to a png image
	t.render(out_file, w=width, units="mm", tree_style=ts)

    #Export the xml file
	project = Phyloxml()
	phylo = phyloxml.PhyloxmlTree(newick=t.write(format=0, features=[]))
	phylo.phyloxml_phylogeny.set_name(title)
	project.add_phylogeny(phylo)
	project.export(open(out_file_xml,'w'))
	
	#Show the tree interactively
	#t.show(tree_style=ts)


if __name__ == "__main__":
	main(sys.argv[1:])



