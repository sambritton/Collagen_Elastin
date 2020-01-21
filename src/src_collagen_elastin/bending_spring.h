#ifndef TORSION_FORCE_H_
#define TORSION_FORCE_H_

struct NodeInfoVecs;
struct BendInfoVecs;
struct GeneralParams;

void calc_bending_spring_force(
	NodeInfoVecs& nodeInfoVecs,
	BendInfoVecs& bendInfoVecs,
	GeneralParams& generalParams);



#endif
