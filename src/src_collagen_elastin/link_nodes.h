#ifndef LINK_NODES_H_
#define LINK_NODES_H_

struct NodeInfoVecs;
struct EdgeInfoVecs;
struct AuxVecs;
struct GeneralParams;

void link_nodes(
	NodeInfoVecs& nodeInfoVecs,
	EdgeInfoVecs& edgeInfoVecs,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);



#endif
