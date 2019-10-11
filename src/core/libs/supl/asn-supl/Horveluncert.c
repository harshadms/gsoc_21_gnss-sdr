/*
 * Generated by asn1c-0.9.22 (http://lionet.info/asn1c)
 * From ASN.1 module "ULP-Components"
 *     found in "../supl-common.asn"
 */

#include "Horveluncert.h"

static int memb_bearing_constraint_1(asn_TYPE_descriptor_t *td,
    const void *sptr,
    asn_app_constraint_failed_f *ctfailcb,
    void *app_key)
{
    const BIT_STRING_t *st = (const BIT_STRING_t *)sptr;
    size_t size;

    if (!sptr)
        {
            _ASN_CTFAIL(app_key, td, sptr, "%s: value not given (%s:%d)",
                td->name, __FILE__, __LINE__);
            return -1;
        }

    if (st->size > 0)
        {
            /* Size in bits */
            size = 8 * st->size - (st->bits_unused & 0x07);
        }
    else
        {
            size = 0;
        }

    if (size == 9)
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

static int memb_horspeed_constraint_1(asn_TYPE_descriptor_t *td,
    const void *sptr,
    asn_app_constraint_failed_f *ctfailcb,
    void *app_key)
{
    const BIT_STRING_t *st = (const BIT_STRING_t *)sptr;
    size_t size;

    if (!sptr)
        {
            _ASN_CTFAIL(app_key, td, sptr, "%s: value not given (%s:%d)",
                td->name, __FILE__, __LINE__);
            return -1;
        }

    if (st->size > 0)
        {
            /* Size in bits */
            size = 8 * st->size - (st->bits_unused & 0x07);
        }
    else
        {
            size = 0;
        }

    if (size == 16)
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

static int memb_uncertspeed_constraint_1(asn_TYPE_descriptor_t *td,
    const void *sptr,
    asn_app_constraint_failed_f *ctfailcb,
    void *app_key)
{
    const BIT_STRING_t *st = (const BIT_STRING_t *)sptr;
    size_t size;

    if (!sptr)
        {
            _ASN_CTFAIL(app_key, td, sptr, "%s: value not given (%s:%d)",
                td->name, __FILE__, __LINE__);
            return -1;
        }

    if (st->size > 0)
        {
            /* Size in bits */
            size = 8 * st->size - (st->bits_unused & 0x07);
        }
    else
        {
            size = 0;
        }

    if (size == 8)
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

static asn_per_constraints_t ASN_PER_MEMB_BEARING_CONSTR_2 = {
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    {APC_CONSTRAINED, 0, 0, 9, 9} /* (SIZE(9..9)) */,
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_HORSPEED_CONSTR_3 = {
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    {APC_CONSTRAINED, 0, 0, 16, 16} /* (SIZE(16..16)) */,
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_UNCERTSPEED_CONSTR_4 = {
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    {APC_CONSTRAINED, 0, 0, 8, 8} /* (SIZE(8..8)) */,
    0,
    0 /* No PER value map */
};
static asn_TYPE_member_t asn_MBR_Horveluncert_1[] = {
    {ATF_NOFLAGS, 0, offsetof(struct Horveluncert, bearing),
        (ASN_TAG_CLASS_CONTEXT | (0 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_BIT_STRING, memb_bearing_constraint_1,
        &ASN_PER_MEMB_BEARING_CONSTR_2, 0, "bearing"},
    {ATF_NOFLAGS, 0, offsetof(struct Horveluncert, horspeed),
        (ASN_TAG_CLASS_CONTEXT | (1 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_BIT_STRING, memb_horspeed_constraint_1,
        &ASN_PER_MEMB_HORSPEED_CONSTR_3, 0, "horspeed"},
    {ATF_NOFLAGS, 0, offsetof(struct Horveluncert, uncertspeed),
        (ASN_TAG_CLASS_CONTEXT | (2 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_BIT_STRING, memb_uncertspeed_constraint_1,
        &ASN_PER_MEMB_UNCERTSPEED_CONSTR_4, 0, "uncertspeed"},
};
static ber_tlv_tag_t asn_DEF_Horveluncert_tags_1[] = {
    (ASN_TAG_CLASS_UNIVERSAL | (16 << 2))};
static asn_TYPE_tag2member_t asn_MAP_Horveluncert_tag2el_1[] = {
    {(ASN_TAG_CLASS_CONTEXT | (0 << 2)), 0, 0, 0}, /* bearing at 245 */
    {(ASN_TAG_CLASS_CONTEXT | (1 << 2)), 1, 0, 0}, /* horspeed at 246 */
    {(ASN_TAG_CLASS_CONTEXT | (2 << 2)), 2, 0, 0}  /* uncertspeed at 247 */
};
static asn_SEQUENCE_specifics_t asn_SPC_Horveluncert_specs_1 = {
    sizeof(struct Horveluncert),
    offsetof(struct Horveluncert, _asn_ctx),
    asn_MAP_Horveluncert_tag2el_1,
    3, /* Count of tags in the map */
    0,
    0,
    0, /* Optional elements (not needed) */
    2, /* Start extensions */
    4  /* Stop extensions */
};
asn_TYPE_descriptor_t asn_DEF_Horveluncert = {
    "Horveluncert",
    "Horveluncert",
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
    asn_DEF_Horveluncert_tags_1,
    sizeof(asn_DEF_Horveluncert_tags_1) /
        sizeof(asn_DEF_Horveluncert_tags_1[0]), /* 1 */
    asn_DEF_Horveluncert_tags_1,                /* Same as above */
    sizeof(asn_DEF_Horveluncert_tags_1) /
        sizeof(asn_DEF_Horveluncert_tags_1[0]), /* 1 */
    0,                                          /* No PER visible constraints */
    asn_MBR_Horveluncert_1,
    3,                            /* Elements count */
    &asn_SPC_Horveluncert_specs_1 /* Additional specs */
};
