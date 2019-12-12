#ifndef PARAMS_CALC_H_
#define PARAMS_CALC_H_


struct NodeInfoVecs;
struct PltInfoVecs;
struct DomainParams;
struct GeneralParams;

void Params_Calc(
    WLCInfoVecs& wlcInfoVecs,
    NodeInfoVecs& nodeInfoVecs,
    GeneralParams& generalParams,
    PltInfoVecs& pltInfoVecs);


								
#endif 