#!/usr/bin/python

#reference tree plot
import networkx as nx
import matplotlib.pyplot as plt
import re
import scipy

try:
    from networkx import graphviz_layout
except ImportError:
    raise ImportError("This example needs Graphviz and either PyGraphviz or Pydot")

def matching(data, search):
	matches = []
	for author, keyword, year in data:
		if author == search:
			matches.append(keyword)
	return matches

######importing bibiography data
from pybtex.database.input import bibtex

parser = bibtex.Parser()
bib_data = parser.parse_file('/home/jeldridge/test_lib.bib')

dataset = []
paired_values = []
for value in bib_data.entries:
	value_author = []
	author_set = re.split(' and ',bib_data.entries[value].fields['author'])
	for things in author_set:
		comma = re.split(',',things)
		a = len(comma[0])+3
		name_and_initial = things[:a]
		value_author.append(name_and_initial)
		paired_values.append([name_and_initial,value,bib_data.entries[value].fields['year']])
	dataset.append([value,value_author,bib_data.entries[value].fields['year']])

authors = []
for item in dataset:
	complete = item[1]  
	for value in complete:
		authors.append(value)
	
#most prolific authors in the field - if the damn script worked!
#plt.hist(authors, bins=len(authors), align="left", histtype="step", color="black")	

#print dataset
#has format: [keyvalue, authors, year] 
node_years = []

for value in dataset:
	node_years.append([value[0],int(value[2])])
	
edges = []
for author, keyword, year in paired_values:
	pairs = matching(paired_values, author)
	for element in pairs:
		edges.append([keyword,element])

nodes_years = sorted(node_years, key=lambda entry: entry[1])

node_list = []
year_list = []

##########turing each key_value into a node
G=nx.Graph()
for nodes, years in node_years:
	G.add_node(nodes, date=years)
	node_list.append(nodes)
	year_list.append(years)
G.add_edges_from(edges)
#print(G.nodes())
#print(G.edges())

####scaling lines to year
single_years = sorted(list(set(year_list)))

##positioning options
#pos=nx.graphviz_layout(G,prog='twopi',args='')
#pos=nx.graphviz_layout(G,prog="neato", overlap="scale")
pos=nx.spring_layout(G)
#pos=nx.random_layout(G)

for n in pos:
	x,y = pos[n]
	y = (single_years.index(G.node[n]['date']))
	#y = years[a]#-min(years))/(max(years)-min(years))
	pos[n] = x,y
	
plt.figure(figsize=(12,12))
#plt.axis('equal')
nx.draw(G,pos,node_size=50,alpha=0.5,node_color="blue", font_size=8)

#plt.savefig('references_map.png')
plt.show()	


	




