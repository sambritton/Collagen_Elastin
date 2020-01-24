#ifndef PARAMS_CALC_H_
#define PARAMS_CALC_H_


struct NodeInfoVecs;
struct PltInfoVecs;
struct DomainParams;
struct GeneralParams;

void params_calc(
    EdgeInfoVecs& edgeInfoVecs,
    NodeInfoVecs& nodeInfoVecs,
    GeneralParams& generalParams,
    PltInfoVecs& pltInfoVecs);


								
#endif 