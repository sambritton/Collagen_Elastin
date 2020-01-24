#ifndef ADVANCE_POSITIONS_H_
#define ADVANCE_POSITIONS_H_

struct NodeInfoVecs;
struct GeneralParams;
struct EdgeInfoVecs;
struct RandVecs;

void advance_positions(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	EdgeInfoVecs& edgeInfoVecs,
	RandVecs& randVecs);



#endif
