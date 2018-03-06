/*
 * Generated by asn1c-0.9.29 (http://lionet.info/asn1c)
 * From ASN.1 module "SUPL-TRIGGERED-START"
 * 	found in "../ulp.asn1"
 * 	`asn1c -S ../../skeletons -pdu=ULP-PDU -pdu=SUPLINIT -fcompound-names -no-gen-OER`
 */

#ifndef _AreaIdList_H_
#define _AreaIdList_H_


#include "asn_application.h"

/* Including external dependencies */
#include "AreaIdSet.h"
#include "AreaIdSetType.h"
#include "constr_SEQUENCE.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* Forward declarations */
    struct GeoAreaMappingList;

    /* AreaIdList */
    typedef struct AreaIdList
    {
        AreaIdSet_t areaIdSet;
        AreaIdSetType_t *areaIdSetType /* OPTIONAL */;
        struct GeoAreaMappingList *geoAreaMappingList /* OPTIONAL */;

        /* Context for parsing across buffer boundaries */
        asn_struct_ctx_t _asn_ctx;
    } AreaIdList_t;

    /* Implementation */
    extern asn_TYPE_descriptor_t asn_DEF_AreaIdList;
    extern asn_SEQUENCE_specifics_t asn_SPC_AreaIdList_specs_1;
    extern asn_TYPE_member_t asn_MBR_AreaIdList_1[3];

#ifdef __cplusplus
}
#endif

/* Referred external types */
#include "GeoAreaMappingList.h"

#endif /* _AreaIdList_H_ */
#include "asn_internal.h"