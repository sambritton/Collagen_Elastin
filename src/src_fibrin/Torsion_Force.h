#ifndef TORSION_FORCE_H_
#define TORSION_FORCE_H_

struct NodeInfoVecs;
struct TorsionInfoVecs;
struct GeneralParams;

void Torsion_Force(
	NodeInfoVecs& nodeInfoVecs,
	TorsionInfoVecs& torsionInfoVecs,
	GeneralParams& generalParams);


								
#endif 