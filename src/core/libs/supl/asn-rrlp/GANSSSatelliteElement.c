/*
 * Generated by asn1c-0.9.22 (http://lionet.info/asn1c)
 * From ASN.1 module "RRLP-Components"
 *     found in "../rrlp-components.asn"
 */

#include "GANSSSatelliteElement.h"

static int memb_svHealth_constraint_1(asn_TYPE_descriptor_t *td,
    const void *sptr,
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

    if ((value >= -7 && value <= 13))
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

static int memb_iod_constraint_1(asn_TYPE_descriptor_t *td, const void *sptr,
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

    if ((value >= 0 && value <= 1023))
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

static asn_per_constraints_t ASN_PER_MEMB_SV_HEALTH_CONSTR_3 = {
    {APC_CONSTRAINED, 5, 5, -7, 13} /* (-7..13) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_IOD_CONSTR_4 = {
    {APC_CONSTRAINED, 10, 10, 0, 1023} /* (0..1023) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_TYPE_member_t asn_MBR_GANSSSatelliteElement_1[] = {
    {ATF_NOFLAGS, 0, offsetof(struct GANSSSatelliteElement, svID),
        (ASN_TAG_CLASS_CONTEXT | (0 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_SVID, 0,                       /* Defer constraints checking to the member type */
        0,                                      /* No PER visible constraints */
        0, "svID"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSSatelliteElement, svHealth),
        (ASN_TAG_CLASS_CONTEXT | (1 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_svHealth_constraint_1,
        &ASN_PER_MEMB_SV_HEALTH_CONSTR_3, 0, "svHealth"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSSatelliteElement, iod),
        (ASN_TAG_CLASS_CONTEXT | (2 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_iod_constraint_1, &ASN_PER_MEMB_IOD_CONSTR_4,
        0, "iod"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSSatelliteElement, ganssClockModel),
        (ASN_TAG_CLASS_CONTEXT | (3 << 2)), +1, /* EXPLICIT tag at current level */
        &asn_DEF_GANSSClockModel,
        0, /* Defer constraints checking to the member type */
        0, /* No PER visible constraints */
        0, "ganssClockModel"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSSatelliteElement, ganssOrbitModel),
        (ASN_TAG_CLASS_CONTEXT | (4 << 2)), +1, /* EXPLICIT tag at current level */
        &asn_DEF_GANSSOrbitModel,
        0, /* Defer constraints checking to the member type */
        0, /* No PER visible constraints */
        0, "ganssOrbitModel"},
};
static ber_tlv_tag_t asn_DEF_GANSSSatelliteElement_tags_1[] = {
    (ASN_TAG_CLASS_UNIVERSAL | (16 << 2))};
static asn_TYPE_tag2member_t asn_MAP_GANSSSatelliteElement_tag2el_1[] = {
    {(ASN_TAG_CLASS_CONTEXT | (0 << 2)), 0, 0, 0}, /* svID at 1238 */
    {(ASN_TAG_CLASS_CONTEXT | (1 << 2)), 1, 0, 0}, /* svHealth at 1239 */
    {(ASN_TAG_CLASS_CONTEXT | (2 << 2)), 2, 0, 0}, /* iod at 1240 */
    {(ASN_TAG_CLASS_CONTEXT | (3 << 2)), 3, 0, 0}, /* ganssClockModel at 1241 */
    {(ASN_TAG_CLASS_CONTEXT | (4 << 2)), 4, 0, 0}  /* ganssOrbitModel at 1242 */
};
static asn_SEQUENCE_specifics_t asn_SPC_GANSSSatelliteElement_specs_1 = {
    sizeof(struct GANSSSatelliteElement),
    offsetof(struct GANSSSatelliteElement, _asn_ctx),
    asn_MAP_GANSSSatelliteElement_tag2el_1,
    5, /* Count of tags in the map */
    0,
    0,
    0, /* Optional elements (not needed) */
    4, /* Start extensions */
    6  /* Stop extensions */
};
asn_TYPE_descriptor_t asn_DEF_GANSSSatelliteElement = {
    "GANSSSatelliteElement",
    "GANSSSatelliteElement",
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
    asn_DEF_GANSSSatelliteElement_tags_1,
    sizeof(asn_DEF_GANSSSatelliteElement_tags_1) /
        sizeof(asn_DEF_GANSSSatelliteElement_tags_1[0]), /* 1 */
    asn_DEF_GANSSSatelliteElement_tags_1,                /* Same as above */
    sizeof(asn_DEF_GANSSSatelliteElement_tags_1) /
        sizeof(asn_DEF_GANSSSatelliteElement_tags_1[0]), /* 1 */
    0,                                                   /* No PER visible constraints */
    asn_MBR_GANSSSatelliteElement_1,
    5,                                     /* Elements count */
    &asn_SPC_GANSSSatelliteElement_specs_1 /* Additional specs */
};
