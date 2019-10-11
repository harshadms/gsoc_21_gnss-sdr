/*
 * Generated by asn1c-0.9.22 (http://lionet.info/asn1c)
 * From ASN.1 module "RRLP-Components"
 *     found in "../rrlp-components.asn"
 */

#include "GANSSUTCModel.h"

static int memb_ganssUtcA1_constraint_1(asn_TYPE_descriptor_t *td,
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

    if ((value >= -8388608 && value <= 8388607))
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

static int memb_ganssUtcA0_constraint_1(asn_TYPE_descriptor_t *td,
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

    if ((value >= (-2147483647L - 1) && value <= 2147483647))
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

static int memb_ganssUtcTot_constraint_1(asn_TYPE_descriptor_t *td,
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

    if ((value >= 0 && value <= 255))
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

static int memb_ganssUtcWNt_constraint_1(asn_TYPE_descriptor_t *td,
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

    if ((value >= 0 && value <= 255))
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

static int memb_ganssUtcDeltaTls_constraint_1(
    asn_TYPE_descriptor_t *td, const void *sptr,
    asn_app_constraint_failed_f *ctfailcb, void *app_key)
{
    long value;

    if (!sptr)
        {
            _ASN_CTFAIL(app_key, td, sptr, "%s: value not given (%s:%d)",
                td->name, __FILE__, __LINE__);
            return -1;
        }

    value = *(const long *)sptr;

    if ((value >= -128 && value <= 127))
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

static int memb_ganssUtcWNlsf_constraint_1(
    asn_TYPE_descriptor_t *td, const void *sptr,
    asn_app_constraint_failed_f *ctfailcb, void *app_key)
{
    long value;

    if (!sptr)
        {
            _ASN_CTFAIL(app_key, td, sptr, "%s: value not given (%s:%d)",
                td->name, __FILE__, __LINE__);
            return -1;
        }

    value = *(const long *)sptr;

    if ((value >= 0 && value <= 255))
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

static int memb_ganssUtcDN_constraint_1(asn_TYPE_descriptor_t *td,
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

    if ((value >= -128 && value <= 127))
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

static int memb_ganssUtcDeltaTlsf_constraint_1(
    asn_TYPE_descriptor_t *td, const void *sptr,
    asn_app_constraint_failed_f *ctfailcb, void *app_key)
{
    long value;

    if (!sptr)
        {
            _ASN_CTFAIL(app_key, td, sptr, "%s: value not given (%s:%d)",
                td->name, __FILE__, __LINE__);
            return -1;
        }

    value = *(const long *)sptr;

    if ((value >= -128 && value <= 127))
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

static asn_per_constraints_t ASN_PER_MEMB_GANSS_UTC_A1_CONSTR_2 = {
    {APC_CONSTRAINED, 24, -1, -8388608, 8388607} /* (-8388608..8388607) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_GANSS_UTC_A0_CONSTR_3 = {
    {APC_CONSTRAINED, 32, -1, (-2147483647L - 1),
        2147483647} /* (-2147483648..2147483647) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_GANSS_UTC_TOT_CONSTR_4 = {
    {APC_CONSTRAINED, 8, 8, 0, 255} /* (0..255) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_GANSS_UTC_W_NT_CONSTR_5 = {
    {APC_CONSTRAINED, 8, 8, 0, 255} /* (0..255) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_GANSS_UTC_DELTA_TLS_CONSTR_6 = {
    {APC_CONSTRAINED, 8, 8, -128, 127} /* (-128..127) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_GANSS_UTC_W_NLSF_CONSTR_7 = {
    {APC_CONSTRAINED, 8, 8, 0, 255} /* (0..255) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_GANSS_UTC_DN_CONSTR_8 = {
    {APC_CONSTRAINED, 8, 8, -128, 127} /* (-128..127) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_per_constraints_t ASN_PER_MEMB_GANSS_UTC_DELTA_TLSF_CONSTR_9 = {
    {APC_CONSTRAINED, 8, 8, -128, 127} /* (-128..127) */,
    {APC_UNCONSTRAINED, -1, -1, 0, 0},
    0,
    0 /* No PER value map */
};
static asn_TYPE_member_t asn_MBR_GANSSUTCModel_1[] = {
    {ATF_NOFLAGS, 0, offsetof(struct GANSSUTCModel, ganssUtcA1),
        (ASN_TAG_CLASS_CONTEXT | (0 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_ganssUtcA1_constraint_1,
        &ASN_PER_MEMB_GANSS_UTC_A1_CONSTR_2, 0, "ganssUtcA1"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSUTCModel, ganssUtcA0),
        (ASN_TAG_CLASS_CONTEXT | (1 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_ganssUtcA0_constraint_1,
        &ASN_PER_MEMB_GANSS_UTC_A0_CONSTR_3, 0, "ganssUtcA0"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSUTCModel, ganssUtcTot),
        (ASN_TAG_CLASS_CONTEXT | (2 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_ganssUtcTot_constraint_1,
        &ASN_PER_MEMB_GANSS_UTC_TOT_CONSTR_4, 0, "ganssUtcTot"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSUTCModel, ganssUtcWNt),
        (ASN_TAG_CLASS_CONTEXT | (3 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_ganssUtcWNt_constraint_1,
        &ASN_PER_MEMB_GANSS_UTC_W_NT_CONSTR_5, 0, "ganssUtcWNt"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSUTCModel, ganssUtcDeltaTls),
        (ASN_TAG_CLASS_CONTEXT | (4 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_ganssUtcDeltaTls_constraint_1,
        &ASN_PER_MEMB_GANSS_UTC_DELTA_TLS_CONSTR_6, 0, "ganssUtcDeltaTls"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSUTCModel, ganssUtcWNlsf),
        (ASN_TAG_CLASS_CONTEXT | (5 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_ganssUtcWNlsf_constraint_1,
        &ASN_PER_MEMB_GANSS_UTC_W_NLSF_CONSTR_7, 0, "ganssUtcWNlsf"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSUTCModel, ganssUtcDN),
        (ASN_TAG_CLASS_CONTEXT | (6 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_ganssUtcDN_constraint_1,
        &ASN_PER_MEMB_GANSS_UTC_DN_CONSTR_8, 0, "ganssUtcDN"},
    {ATF_NOFLAGS, 0, offsetof(struct GANSSUTCModel, ganssUtcDeltaTlsf),
        (ASN_TAG_CLASS_CONTEXT | (7 << 2)), -1, /* IMPLICIT tag at current level */
        &asn_DEF_NativeInteger, memb_ganssUtcDeltaTlsf_constraint_1,
        &ASN_PER_MEMB_GANSS_UTC_DELTA_TLSF_CONSTR_9, 0, "ganssUtcDeltaTlsf"},
};
static ber_tlv_tag_t asn_DEF_GANSSUTCModel_tags_1[] = {
    (ASN_TAG_CLASS_UNIVERSAL | (16 << 2))};
static asn_TYPE_tag2member_t asn_MAP_GANSSUTCModel_tag2el_1[] = {
    {(ASN_TAG_CLASS_CONTEXT | (0 << 2)), 0, 0, 0}, /* ganssUtcA1 at 1382 */
    {(ASN_TAG_CLASS_CONTEXT | (1 << 2)), 1, 0, 0}, /* ganssUtcA0 at 1383 */
    {(ASN_TAG_CLASS_CONTEXT | (2 << 2)), 2, 0, 0}, /* ganssUtcTot at 1384 */
    {(ASN_TAG_CLASS_CONTEXT | (3 << 2)), 3, 0, 0}, /* ganssUtcWNt at 1385 */
    {(ASN_TAG_CLASS_CONTEXT | (4 << 2)), 4, 0,
        0},                                        /* ganssUtcDeltaTls at 1386 */
    {(ASN_TAG_CLASS_CONTEXT | (5 << 2)), 5, 0, 0}, /* ganssUtcWNlsf at 1387 */
    {(ASN_TAG_CLASS_CONTEXT | (6 << 2)), 6, 0, 0}, /* ganssUtcDN at 1388 */
    {(ASN_TAG_CLASS_CONTEXT | (7 << 2)), 7, 0,
        0} /* ganssUtcDeltaTlsf at 1389 */
};
static asn_SEQUENCE_specifics_t asn_SPC_GANSSUTCModel_specs_1 = {
    sizeof(struct GANSSUTCModel),
    offsetof(struct GANSSUTCModel, _asn_ctx),
    asn_MAP_GANSSUTCModel_tag2el_1,
    8, /* Count of tags in the map */
    0,
    0,
    0,  /* Optional elements (not needed) */
    -1, /* Start extensions */
    -1  /* Stop extensions */
};
asn_TYPE_descriptor_t asn_DEF_GANSSUTCModel = {
    "GANSSUTCModel",
    "GANSSUTCModel",
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
    asn_DEF_GANSSUTCModel_tags_1,
    sizeof(asn_DEF_GANSSUTCModel_tags_1) /
        sizeof(asn_DEF_GANSSUTCModel_tags_1[0]), /* 1 */
    asn_DEF_GANSSUTCModel_tags_1,                /* Same as above */
    sizeof(asn_DEF_GANSSUTCModel_tags_1) /
        sizeof(asn_DEF_GANSSUTCModel_tags_1[0]), /* 1 */
    0,                                           /* No PER visible constraints */
    asn_MBR_GANSSUTCModel_1,
    8,                             /* Elements count */
    &asn_SPC_GANSSUTCModel_specs_1 /* Additional specs */
};
