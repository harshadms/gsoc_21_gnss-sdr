/*
 * Generated by asn1c-0.9.22 (http://lionet.info/asn1c)
 * From ASN.1 module "RRLP-Components"
 *     found in "../rrlp-components.asn"
 */

#ifndef _SgnTypeElement_H_
#define _SgnTypeElement_H_

#include <asn_application.h>

/* Including external dependencies */
#include "GANSSSignalID.h"
#include "SeqOfDGANSSSgnElement.h"
#include <NativeInteger.h>
#include <constr_SEQUENCE.h>

#ifdef __cplusplus
extern "C"
{
#endif

    /* SgnTypeElement */
    typedef struct SgnTypeElement
    {
        GANSSSignalID_t *ganssSignalID /* OPTIONAL */;
        long ganssStatusHealth;
        SeqOfDGANSSSgnElement_t dganssSgnList;

        /* Context for parsing across buffer boundaries */
        asn_struct_ctx_t _asn_ctx;
    } SgnTypeElement_t;

    /* Implementation */
    extern asn_TYPE_descriptor_t asn_DEF_SgnTypeElement;

#ifdef __cplusplus
}
#endif

#endif /* _SgnTypeElement_H_ */
#include <asn_internal.h>
