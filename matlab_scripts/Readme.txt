The set of functions generates a mesh with given degree, alignment and length distributions for a given size.
Usage: 
generate_graph.m | builds xml file and Graph object in two steps.
1. generate_initial_structure.m | builds initial set of edges matching degree distributions
2. length_alignment_matcher.m | alters point locations to match length and alignment distributions

generate_graph_parameters.m | takes a Graph object and builds the xml file for the c++ code. 

Use as follows: "generate_graph_parameters(generate_graph())"