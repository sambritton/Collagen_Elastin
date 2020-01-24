#ifndef BENDING_SPRING_H_
#define BENDING_SPRING_H_

struct NodeInfoVecs;
struct BendInfoVecs;
struct GeneralParams;

void calc_bending_spring_force(
	NodeInfoVecs& nodeInfoVecs,
	BendInfoVecs& bendInfoVecs,
	GeneralParams& generalParams);



#endif
