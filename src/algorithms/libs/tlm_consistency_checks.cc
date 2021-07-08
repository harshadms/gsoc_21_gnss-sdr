/*!
 * \file tlm_consistency_checks.cc
 * \brief Library of anti-spoofing functions
 * \author Harshad Sathaye sathaye.h(at)northeastern.edu
 *
 * -----------------------------------------------------------------------------
 *
 * GNSS-SDR is a Global Navigation Satellite System software-defined receiver.
 * This file is part of GNSS-SDR.
 *
 * Copyright (C) 2010-2021  (see AUTHORS file for a list of contributors)
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * -----------------------------------------------------------------------------
 */

#include "spoofing_detector.h"
#include <cmath>  // for floor, fmod, rint, ceil

TLMConsistencyChecks::TLMConsistencyChecks()
{
}

TLMConsistencyChecks::TLMConsistencyChecks(const TLMConsistencyChecksConf *conf_)
{
    d_check_RX_clock = conf_->check_RX_clock;
    d_check_TOW = conf_->check_TOW;
    DLOG(INFO) << "TLM_CHECKS: RX_clock" << d_check_RX_clock;
    DLOG(INFO) << "TLM_CHECKS: TOW" << d_check_TOW;

    d_first_record = true;
}

void TLMConsistencyChecks::check_RX_clock()
{
    if (d_first_record)
        {
            set_old_clock();
            set_lkg_clock(false);
            d_first_record = false;
            return;
        }

    check_clock_jump();
}

void TLMConsistencyChecks::check_clock_jump()
{
    double sample_diff = d_new_clock.sample_counter - d_old_clock.sample_counter;
    double sample_diff_time = sample_diff * 1000 / d_gnss_synchro.fs;

    uint32_t tow_diff = d_new_clock.tow - d_old_clock.tow;

    uint32_t time_diff = abs(tow_diff - sample_diff_time);

    if (time_diff > 20)
        {
            DLOG(INFO) << "TLM_CHECKS:  PRN " << PRN << " Sample diff - " << sample_diff_time << " TOW diff - " << tow_diff;
        }

    DLOG(INFO) << "TLM_CHECKS:  PRN " << PRN << " Sample diff - " << sample_diff_time << " TOW diff - " << tow_diff;
    set_old_clock();
    set_lkg_clock(false);
}

void TLMConsistencyChecks::update_clock_info(uint64_t sample_counter, uint32_t tow, uint32_t wn)
{
    d_new_clock.sample_counter = sample_counter;
    d_new_clock.tow = tow;
    d_new_clock.wn = wn;
    //DLOG(INFO) << "TLM_CHECKS: Received new clock info - Channel " << channel_id << " : PRN " << PRN << " = " << sample_counter << " " << tow;
}

void TLMConsistencyChecks::set_old_clock()
{
    d_old_clock.sample_counter = d_new_clock.sample_counter;
    d_old_clock.tow = d_new_clock.tow;
    d_old_clock.wn = d_new_clock.wn;
}

void TLMConsistencyChecks::set_lkg_clock(bool set_old)
{
    if (set_old)
        {
            d_lkg_clock.sample_counter = d_old_clock.sample_counter;
            d_lkg_clock.tow = d_old_clock.tow;
            d_lkg_clock.wn = d_old_clock.wn;
        }
    else
        {
            d_lkg_clock.sample_counter = d_new_clock.sample_counter;
            d_lkg_clock.tow = d_new_clock.tow;
            d_lkg_clock.wn = d_new_clock.wn;
        }

    // Optional DLOG message
}