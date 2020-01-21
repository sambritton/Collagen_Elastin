#ifndef BUCKET_SCHEME_H_
#define BUCKET_SCHEME_H_


struct NodeInfoVecs;
struct GeneralParams;
struct DomainParams;
struct EdgeInfoVecs;
struct BendInfoVecs;
struct AuxVecs;

void init_dim_general(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);

void init_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);

void build_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);

void extend_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);

#endif
