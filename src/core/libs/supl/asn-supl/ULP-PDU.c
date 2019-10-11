/*
 * Generated by asn1c-0.9.22 (http://lionet.info/asn1c)
 * From ASN.1 module "ULP"
 *     found in "../supl-ulp.asn"
 */

#include "ULP-PDU.h"

static int memb_length_constraint_1(asn_TYPE_descriptor_t *td, const void *sptr,
    asn_app_constraint_failed_f *ctfailcb,
    void *app_key)
{
    long value;

    if (!sptr)
        {
            _ASN_CTFAIL(app_key, td, sptr, "%s: value not given (%s:%d)",
                td->name, __FILE__, __LINE__);
            return -1;
        }

    value = *(const long *)sptr;

    if ((value >= 0 && value <= 65535))
        {
            /* Constraint check succeeded */
            return 0;
        }
    else
        {
            _ASN_CTFAIL(app_key, td, sptr, "%s: constraint failed (%s:%d)",
                td->name, __FILE__, __LINE__);
            return -1;
        }
}

static asn_per_constraints_t ASN_PER_MEMB_LENGTH_CONSTR_2 = {
    {APC_CONSTRAINED, 16, 16, 0, 65535} /* (0..65535) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_TYPE_member_t asn_MBR_ULP_PDU_1[] = {
    {ATF_NOFLAGS, 0, offsetof(struct ULP_PDU, length),
        (ASN_TAG_CLASS_CONTEXT | (0 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_length_constraint_1,
        &ASN_PER_MEMB_LENGTH_CONSTR_2, 0, "length"},
    {ATF_NOFLAGS, 0, offsetof(struct ULP_PDU, version),
        (ASN_TAG_CLASS_CONTEXT | (1 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_Version, 0,                    /* Defer constraints checking to the member type */
        0,                                      /* No PER visible constraints */
        0, "version"},
    {ATF_NOFLAGS, 0, offsetof(struct ULP_PDU, sessionID),
        (ASN_TAG_CLASS_CONTEXT | (2 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_SessionID, 0,                  /* Defer constraints checking to the member type */
        0,                                      /* No PER visible constraints */
        0, "sessionID"},
    {ATF_NOFLAGS, 0, offsetof(struct ULP_PDU, message),
        (ASN_TAG_CLASS_CONTEXT | (3 << 2)), +1, /* EXPLICIT tag at current level */
        &asn_DEF_UlpMessage, 0,                 /* Defer constraints checking to the member type */
        0,                                      /* No PER visible constraints */
        0, "message"},
};
static ber_tlv_tag_t asn_DEF_ULP_PDU_tags_1[] = {
    (ASN_TAG_CLASS_UNIVERSAL | (16 << 2))};
static asn_TYPE_tag2member_t asn_MAP_ULP_PDU_tag2el_1[] = {
    {(ASN_TAG_CLASS_CONTEXT | (0 << 2)), 0, 0, 0}, /* length at 31 */
    {(ASN_TAG_CLASS_CONTEXT | (1 << 2)), 1, 0, 0}, /* version at 32 */
    {(ASN_TAG_CLASS_CONTEXT | (2 << 2)), 2, 0, 0}, /* sessionID at 33 */
    {(ASN_TAG_CLASS_CONTEXT | (3 << 2)), 3, 0, 0}  /* message at 34 */
};
static asn_SEQUENCE_specifics_t asn_SPC_ULP_PDU_specs_1 = {
    sizeof(struct ULP_PDU),
    offsetof(struct ULP_PDU, _asn_ctx),
    asn_MAP_ULP_PDU_tag2el_1,
    4, /* Count of tags in the map */
    0,
    0,
    0,  /* Optional elements (not needed) */
    -1, /* Start extensions */
    -1  /* Stop extensions */
};
asn_TYPE_descriptor_t asn_DEF_ULP_PDU = {
    "ULP-PDU",
    "ULP-PDU",
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
    asn_DEF_ULP_PDU_tags_1,
    sizeof(asn_DEF_ULP_PDU_tags_1) / sizeof(asn_DEF_ULP_PDU_tags_1[0]), /* 1 */
    asn_DEF_ULP_PDU_tags_1,                                             /* Same as above */
    sizeof(asn_DEF_ULP_PDU_tags_1) / sizeof(asn_DEF_ULP_PDU_tags_1[0]), /* 1 */
    0,                                                                  /* No PER visible constraints */
    asn_MBR_ULP_PDU_1,
    4,                       /* Elements count */
    &asn_SPC_ULP_PDU_specs_1 /* Additional specs */
};
