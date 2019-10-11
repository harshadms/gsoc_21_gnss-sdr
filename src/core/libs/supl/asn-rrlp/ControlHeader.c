/*
 * Generated by asn1c-0.9.22 (http://lionet.info/asn1c)
 * From ASN.1 module "RRLP-Components"
 *     found in "../rrlp-components.asn"
 */

#include "ControlHeader.h"

static asn_TYPE_member_t asn_MBR_ControlHeader_1[] = {
    {ATF_POINTER, 9, offsetof(struct ControlHeader, referenceTime),
        (ASN_TAG_CLASS_CONTEXT | (0 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_ReferenceTime,
        0, /* Defer constraints checking to the member type */
        0, /* No PER visible constraints */
        0, "referenceTime"},
    {ATF_POINTER, 8, offsetof(struct ControlHeader, refLocation),
        (ASN_TAG_CLASS_CONTEXT | (1 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_RefLocation,
        0, /* Defer constraints checking to the member type */
        0, /* No PER visible constraints */
        0, "refLocation"},
    {ATF_POINTER, 7, offsetof(struct ControlHeader, dgpsCorrections),
        (ASN_TAG_CLASS_CONTEXT | (2 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_DGPSCorrections,
        0, /* Defer constraints checking to the member type */
        0, /* No PER visible constraints */
        0, "dgpsCorrections"},
    {ATF_POINTER, 6, offsetof(struct ControlHeader, navigationModel),
        (ASN_TAG_CLASS_CONTEXT | (3 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NavigationModel,
        0, /* Defer constraints checking to the member type */
        0, /* No PER visible constraints */
        0, "navigationModel"},
    {ATF_POINTER, 5, offsetof(struct ControlHeader, ionosphericModel),
        (ASN_TAG_CLASS_CONTEXT | (4 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_IonosphericModel,
        0, /* Defer constraints checking to the member type */
        0, /* No PER visible constraints */
        0, "ionosphericModel"},
    {ATF_POINTER, 4, offsetof(struct ControlHeader, utcModel),
        (ASN_TAG_CLASS_CONTEXT | (5 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_UTCModel, 0,                   /* Defer constraints checking to the member type */
        0,                                      /* No PER visible constraints */
        0, "utcModel"},
    {ATF_POINTER, 3, offsetof(struct ControlHeader, almanac),
        (ASN_TAG_CLASS_CONTEXT | (6 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_Almanac, 0,                    /* Defer constraints checking to the member type */
        0,                                      /* No PER visible constraints */
        0, "almanac"},
    {ATF_POINTER, 2, offsetof(struct ControlHeader, acquisAssist),
        (ASN_TAG_CLASS_CONTEXT | (7 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_AcquisAssist,
        0, /* Defer constraints checking to the member type */
        0, /* No PER visible constraints */
        0, "acquisAssist"},
    {ATF_POINTER, 1, offsetof(struct ControlHeader, realTimeIntegrity),
        (ASN_TAG_CLASS_CONTEXT | (8 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_SeqOf_BadSatelliteSet,
        0, /* Defer constraints checking to the member type */
        0, /* No PER visible constraints */
        0, "realTimeIntegrity"},
};
static int asn_MAP_ControlHeader_oms_1[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
static ber_tlv_tag_t asn_DEF_ControlHeader_tags_1[] = {
    (ASN_TAG_CLASS_UNIVERSAL | (16 << 2))};
static asn_TYPE_tag2member_t asn_MAP_ControlHeader_tag2el_1[] = {
    {(ASN_TAG_CLASS_CONTEXT | (0 << 2)), 0, 0, 0}, /* referenceTime at 574 */
    {(ASN_TAG_CLASS_CONTEXT | (1 << 2)), 1, 0, 0}, /* refLocation at 575 */
    {(ASN_TAG_CLASS_CONTEXT | (2 << 2)), 2, 0, 0}, /* dgpsCorrections at 576 */
    {(ASN_TAG_CLASS_CONTEXT | (3 << 2)), 3, 0, 0}, /* navigationModel at 577 */
    {(ASN_TAG_CLASS_CONTEXT | (4 << 2)), 4, 0, 0}, /* ionosphericModel at 578 */
    {(ASN_TAG_CLASS_CONTEXT | (5 << 2)), 5, 0, 0}, /* utcModel at 579 */
    {(ASN_TAG_CLASS_CONTEXT | (6 << 2)), 6, 0, 0}, /* almanac at 580 */
    {(ASN_TAG_CLASS_CONTEXT | (7 << 2)), 7, 0, 0}, /* acquisAssist at 581 */
    {(ASN_TAG_CLASS_CONTEXT | (8 << 2)), 8, 0, 0}  /* realTimeIntegrity at 582 */
};
static asn_SEQUENCE_specifics_t asn_SPC_ControlHeader_specs_1 = {
    sizeof(struct ControlHeader),
    offsetof(struct ControlHeader, _asn_ctx),
    asn_MAP_ControlHeader_tag2el_1,
    9,                           /* Count of tags in the map */
    asn_MAP_ControlHeader_oms_1, /* Optional members */
    9,
    0,  /* Root/Additions */
    -1, /* Start extensions */
    -1  /* Stop extensions */
};
asn_TYPE_descriptor_t asn_DEF_ControlHeader = {
    "ControlHeader",
    "ControlHeader",
    SEQUENCE_free,
    SEQUENCE_print,
    SEQUENCE_constraint,
    SEQUENCE_decode_ber,
    SEQUENCE_encode_der,
    SEQUENCE_decode_xer,
    SEQUENCE_encode_xer,
    SEQUENCE_decode_uper,
    SEQUENCE_encode_uper,
    0, /* Use generic outmost tag fetcher */
    asn_DEF_ControlHeader_tags_1,
    sizeof(asn_DEF_ControlHeader_tags_1) /
        sizeof(asn_DEF_ControlHeader_tags_1[0]), /* 1 */
    asn_DEF_ControlHeader_tags_1,                /* Same as above */
    sizeof(asn_DEF_ControlHeader_tags_1) /
        sizeof(asn_DEF_ControlHeader_tags_1[0]), /* 1 */
    0,                                           /* No PER visible constraints */
    asn_MBR_ControlHeader_1,
    9,                             /* Elements count */
    &asn_SPC_ControlHeader_specs_1 /* Additional specs */
};
