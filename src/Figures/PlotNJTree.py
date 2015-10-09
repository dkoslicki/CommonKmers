import numpy as np
import h5py
from Bio import Phylo
from Bio.Phylo import NewickIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from ete2 import Tree, faces, TreeStyle, COLOR_SCHEMES, TextFace, BarChartFace, CircleFace, AttrFace, Phyloxml, phyloxml
import math
import sys, getopt
import os

#python PlotNJTree.py -i ../Data/test2-reads.fa-x.txt -D ../Data/ -l -n -o ../Output/test.png -x ../Output/test.png

def main(argv):
	input_file=''
	title='Title'
	label_internal_nodes = False
	label_leaves = False
	out_file=''
	width=750
	out_file_xml=''
	plot_rectangular = False
	common_kmer_data_path=''
	taxonomic_names_on_leaves = False
	try:
		opts, args = getopt.getopt(argv,"h:i:lnrto:w:x:D:",["Help=","InputCommonKmerXFile=","LabelLeaves=", "LabelInternalNodes=","Rectangular=", "TaxonomicNamesOnLeaves=", "OutFile=","Width=","OutFileXML=","CommonKmerDataPath="])
	except getopt.GetoptError:
		print 'Unknown option, call using: ./PlotNJTree.py -i <InputCommonKmerXFile> -D <CommonKmerDataPath> -l <LabelLeavesFlag> -n <LabelInternalNodesFlag> -r <RectangularPlotFlag> -t <TaxonomicNamesOnLeavesFlag> -o <OutFile.png> -x <Outfile.xml> -w <Width>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './PlotNJTree.py -i <InputCommonKmerXFile> -D <CommonKmerDataPath> -l <LabelLeavesFlag> -n <LabelInternalNodesFlag> -r <RectangularPlotFlag> -t <TaxonomicNamesOnLeavesFlag> -o <OutFile.png> -x <Outfile.xml> -w <Width>'
			sys.exit(2)
		elif opt in ("-i", "--InputCommonKmerXFile"):
			input_file = arg
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
		elif opt in ("-D", "--CommonKmerDataPath"):
			common_kmer_data_path = arg
		elif opt in ("-r", "--Rectangular"):
			plot_rectangular = True
		elif opt in ("-t", "--TaxonomicNamesOnLeaves"):
			taxonomic_names_on_leaves = True
	
	
	#Read in the x vector
	fid = open(input_file,'r')
	x = map(lambda y: float(y),fid.readlines())
	fid.close()
	
	#Normalize the x vector
	#x = map(lambda y: y/sum(x),x)
	
	#Read in the taxonomy
	taxonomy = list()
	fid = open(os.path.join(common_kmer_data_path,"Taxonomy.txt"),'r')
	for line in fid:
		taxonomy.append('_'.join(line.split()[0].split("_")[1:])) #Just take the first line of the taxonomy (erasing the taxID)
	fid.close()
	
	#Read in the basis for the ckm matrices
	x_file_names = list()
	fid = open(os.path.join(common_kmer_data_path,"FileNames.txt"),'r')
	for line in fid:
		x_file_names.append(os.path.basename(line.strip()))
	fid.close()
	
	#Read in the common kmer matrix
	f=h5py.File(os.path.join(common_kmer_data_path,'CommonKmerMatrix-30mers.h5'),'r')
	ckm30=np.array(f['common_kmers'],dtype=np.float64)
	f.close()
	f=h5py.File(os.path.join(common_kmer_data_path,'CommonKmerMatrix-50mers.h5'),'r')
	ckm50=np.array(f['common_kmers'],dtype=np.float64)
	f.close()
	ckm30_norm = np.multiply(ckm30,1/np.diag(ckm30))
	ckm50_norm = np.multiply(ckm50,1/np.diag(ckm50))
	num_rows = ckm30_norm.shape[0]
	num_cols = ckm30_norm.shape[1]
	names = x_file_names
	matrix=list()
	for i in range(num_rows):
		matrix.append([.5*(1-.5*ckm30_norm[i,j]-.5*ckm30_norm[j,i])+.5*(1-.5*ckm50_norm[i,j]-.5*ckm50_norm[j,i]) for j in range(i+1)])
	
	#Construct the tree. Note I could use RapidNJ here, but a few tests have shown that the trees that RapidNJ creates are rubbish.
	dm = _DistanceMatrix(names, matrix)
	constructor = DistanceTreeConstructor()
	tree = constructor.nj(dm)
	t=Tree(tree.format('newick'),format=1)
	#tree.format('newick')
	#Phylo.draw_ascii(tree)
	
	#Now I will put internal nodes in a certain phylogenetic distance between the root and a given node.
	#Function to insert a node at a given distance
	def insert_node(t, name_to_insert, insert_above, dist_along):
		insert_at_node = t.search_nodes(name=insert_above)[0]
		parent = (t&insert_above).up
		orig_branch_length = t.get_distance(insert_at_node,parent)
		if orig_branch_length < dist_along:
			raise ValueError("error: dist_along larger than orig_branch_length")
		removed_node = insert_at_node.detach()
		removed_node.dist = orig_branch_length - dist_along
		added_node = parent.add_child(name=name_to_insert, dist=dist_along)
		added_node.add_child(removed_node)
	
	#Function to insert a node some % along a branch
	def insert_hyp_node(t, leaf_name, percent):
		total_dist = t.get_distance(t.name,leaf_name)
		percent_dist = percent*total_dist
		child_node = (t&leaf_name)
		ancestor_node = (t&child_node.name).up
		while t.get_distance(t.name, ancestor_node) > percent_dist:
			child_node = ancestor_node
			ancestor_node = (t&child_node.name).up
		insert_node(t, leaf_name+"_"+str(percent), child_node.name, percent_dist-t.get_distance(t.name, ancestor_node))
	
	#Insert hypothetical nodes
	hyp_node_names = dict()
	cutoffs = [.9,.8,.7,.6,.5,.4,.3,.2,.1]
	cutoffs = map(lambda y: y**1.5,cutoffs)
	for i in range(len(x_file_names)):
		xi = x[i:len(x):len(x_file_names)]
		for j in range(1,len(cutoffs)+1):
			if xi[j]>0:
				insert_hyp_node(t, x_file_names[i], cutoffs[j-1])
				hyp_node_names[x_file_names[i]+"_"+str(cutoffs[j-1])] = [x_file_names[i], cutoffs[j-1], j-1] #in case there are "_" in the file names
				#insert_hyp_node(t, x_file_names[i],.5/t.get_distance(t.name,t&x_file_names[i])*cutoffs[j])
	
	#Now put the bubbles on the nodes
	def layout(node):
		#print(node)
		if node.is_leaf():
			if node.name in x_file_names:
				#make reconstructed bubble
				size = x[x_file_names.index(node.name)]
				F = CircleFace(radius=500*math.sqrt(size), color="RoyalBlue", style="sphere")
				F.border.width = None
				F.opacity = 0.6
				faces.add_face_to_node(F,node, 0, position="branch-right")
				if taxonomic_names_on_leaves:
					nameFace = AttrFace("name", fsize=25, fgcolor='black',text_suffix="_"+taxonomy[x_file_names.index(node.name)])
					faces.add_face_to_node(nameFace, node, 0, position="branch-right")
				else:
					nameFace = AttrFace("name", fsize=25, fgcolor='black')
					faces.add_face_to_node(nameFace, node, 0, position="branch-right")
		elif node.name in hyp_node_names: #Otherwise it's a hypothetical node, just use recon x
			node_base_name = hyp_node_names[node.name][0]
			percent = hyp_node_names[node.name][1]
			if node_base_name in x_file_names:
				idx = hyp_node_names[node.name][2]
				size = x[x_file_names.index(node_base_name)+(idx+1)*len(x_file_names)]
				F = CircleFace(radius=500*math.sqrt(size), color="RoyalBlue", style="sphere")
				F.border.width = None
				F.opacity = 0.6
				faces.add_face_to_node(F,node, 0, position="branch-right")
				#print node
				#print size
			else:
				size=0
		else:
			size=0
		#print(size)
	
	ts = TreeStyle()
	ts.layout_fn = layout
	if plot_rectangular:
		ts.mode = "r"
	else:
		ts.mode = "c"
	ts.show_leaf_name = False
	ts.min_leaf_separation = 50

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
