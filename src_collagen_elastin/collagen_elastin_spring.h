#ifndef COLLAGEN_ELASTIN_SPRING_H_
#define COLLAGEN_ELASTIN_SPRING_H_

struct NodeInfoVecs;
struct EdgeInfoVecs;
struct DomainParams;
struct GeneralParams;

void calc_spring_force(
	NodeInfoVecs& nodeInfoVecs,
	EdgeInfoVecs& edgeInfoVecs,
	GeneralParams& generalParams);

//void get_strain_parameters(
//	NodeInfoVecs& nodeInfoVecs,
//	EdgeInfoVecs& edgeInfoVecs,
//	GeneralParams& generalParams,
//	DomainParams& domainParams);

#endif
