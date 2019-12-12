#ifndef WLC_FORCE_H_
#define WLC_FORCE_H_

struct NodeInfoVecs;
struct WLCInfoVecs;
struct DomainParams;
struct GeneralParams;

void WLC_Force(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams);
	
void GetStrainParameters(NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,  
	GeneralParams& generalParams,
	DomainParams& domainParams);




#endif
