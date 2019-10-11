/*
 * Generated by asn1c-0.9.22 (http://lionet.info/asn1c)
 * From ASN.1 module "RRLP-Components"
 *     found in "../rrlp-components.asn"
 */

#ifndef _MsrPosition_Rsp_H_
#define _MsrPosition_Rsp_H_

#include <asn_application.h>

/* Including external dependencies */
#include "ExtensionContainer.h"
#include <constr_SEQUENCE.h>

#ifdef __cplusplus
extern "C"
{
#endif

    /* Forward declarations */
    struct MultipleSets;
    struct ReferenceIdentity;
    struct OTD_MeasureInfo;
    struct LocationInfo;
    struct GPS_MeasureInfo;
    struct LocationError;
    struct Rel_98_MsrPosition_Rsp_Extension;
    struct Rel_5_MsrPosition_Rsp_Extension;

    /* MsrPosition-Rsp */
    typedef struct MsrPosition_Rsp
    {
        struct MultipleSets *multipleSets /* OPTIONAL */;
        struct ReferenceIdentity *referenceIdentity /* OPTIONAL */;
        struct OTD_MeasureInfo *otd_MeasureInfo /* OPTIONAL */;
        struct LocationInfo *locationInfo /* OPTIONAL */;
        struct GPS_MeasureInfo *gps_MeasureInfo /* OPTIONAL */;
        struct LocationError *locationError /* OPTIONAL */;
        ExtensionContainer_t *extensionContainer /* OPTIONAL */;
        /*
         * This type is extensible,
         * possible extensions are below.
         */
        struct Rel_98_MsrPosition_Rsp_Extension
            *rel_98_MsrPosition_Rsp_Extension /* OPTIONAL */;
        struct Rel_5_MsrPosition_Rsp_Extension
            *rel_5_MsrPosition_Rsp_Extension /* OPTIONAL */;

        /* Context for parsing across buffer boundaries */
        asn_struct_ctx_t _asn_ctx;
    } MsrPosition_Rsp_t;

    /* Implementation */
    extern asn_TYPE_descriptor_t asn_DEF_MsrPosition_Rsp;

#ifdef __cplusplus
}
#endif

/* Referred external types */
#include "GPS-MeasureInfo.h"
#include "LocationError.h"
#include "LocationInfo.h"
#include "MultipleSets.h"
#include "OTD-MeasureInfo.h"
#include "ReferenceIdentity.h"
#include "Rel-5-MsrPosition-Rsp-Extension.h"
#include "Rel-98-MsrPosition-Rsp-Extension.h"

#endif /* _MsrPosition_Rsp_H_ */
#include <asn_internal.h>
