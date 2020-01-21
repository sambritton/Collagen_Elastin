#ifndef BUCKET_SCHEME_H_
#define BUCKET_SCHEME_H_


struct NodeInfoVecs;
struct GeneralParams;
struct ExtensionParams;
struct DomainParams;
struct EdgeInfoVecs;
struct BendInfoVecs;
struct AuxVecs;



void external_force(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	ExtensionParams& extensionParams,
	DomainParams& domainParams);

#endif
