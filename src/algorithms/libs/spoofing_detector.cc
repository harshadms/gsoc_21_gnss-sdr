/*!
 * \file spoofing_detector.cc
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
#include <map>

SpoofingDetector::SpoofingDetector()
{
}

SpoofingDetector::SpoofingDetector(const Pvt_SD_Conf* conf_)
{
    DLOG(INFO) << "Spoofing detector for PVT initialized";
    d_max_jump_distance = conf_->max_jump_distance;
    d_geo_fence_radius = conf_->geo_fence_radius;
    d_velocity_difference = conf_->velocity_difference;
    d_pos_jump_recovery = conf_->pos_jump_recovery;

    d_static_lat = conf_->static_lat;
    d_static_lon = conf_->static_lon;
    d_static_alt = conf_->static_alt;

    DLOG(INFO) << "Static pos check: " << d_static_lat << ", " << d_static_lon;

    d_dump_pos_checks_results = conf_->dump_pos_checks_results;
    d_position_check = conf_->position_check;
    d_static_pos_check = conf_->static_pos_check;

    d_spoofer_score = 0;

    d_old_lat = 0;
    d_old_lon = 0;
    d_old_alt = 0;

    d_last_known_good_lat = 0;
    d_last_known_good_lon = 0;
    d_last_known_good_alt = 0;

    // Used to set old coordinates and last known good location
    d_first_record = true;

    //  Used to decide whether to update LKGL
    d_update_lkgl = true;
}

// ####### Position consistency functions
void SpoofingDetector::check_position_consistency(double lat, double lon, double alt,
    const Gnss_Synchro** in)
{
    // Need to check this for dynamic scenarios.
    if (d_static_pos_check)
        {
            static_pos_check(lat, lon, alt);
        }

    position_jump(lat, lon, alt);

    if (d_first_record)
        {
            d_first_record = false;
        }
}

void SpoofingDetector::compare_velocity()
{
    //DLOG(INFO) << "Compare velocity function";
}

void SpoofingDetector::static_pos_check(double lat, double lon, double alt)
{
    long double distance = calculate_distance(lat, lon, d_static_lat, d_static_lon);
    DLOG(INFO) << "STATIC_POS: static_coords (" << d_static_lat << "," << d_static_lon << ") received (" << lat << "," << lon << ") Distance: " << distance;
    if (distance > d_geo_fence_radius)
        {
            d_score.static_pos_check_score = 1;
            d_spoofer_score = d_score.total_score();
            DLOG(INFO) << "STATIC_POS: Spoofer score: " << d_spoofer_score << " - Distance: " << distance;
        }
    else
        {
            if (d_score.static_pos_check_score == 1)
                {
                    DLOG(INFO) << "STATIC_POS: Received location within geo-fence, resetting static position score";
                }

            d_score.static_pos_check_score = 0;
            d_spoofer_score = d_score.total_score();
        }
}

void SpoofingDetector::position_jump(double new_lat, double new_lon, double new_alt)
{
    double distance_to_lkgl;  // LKGL - Last Known Good Location
    double jump_distance;

    DLOG(INFO) << "POS_JUMP: check";

    if (d_first_record)
        {
            set_old_location(new_lat, new_lon, new_alt);

            // In case of a cold start there is no way to know if the received location is legit. Hence a naive way is to check for spoofer score.
            // The score will be 0 if no other technique detects spoofing
            if (get_spoofer_score() == 0)
                {
                    // Set last known good location to current coordinates
                    set_last_known_good_location(new_lat, new_lon, new_alt);
                    DLOG(INFO) << "POS_JUMP: location update: " << new_lat << ", " << new_lon << ", " << new_alt;
                }

            d_first_record = false;
            return;
        }

    jump_distance = calculate_distance(d_old_lat, d_old_lon, new_lat, new_lon);
    DLOG(INFO) << "POS_JUMP: Old (" << d_old_lat << "," << d_old_lon << ") received (" << new_lat << "," << new_lon << ") Distance: " << jump_distance << " Spoofer score: " << get_spoofer_score();
    if (jump_distance > d_max_jump_distance)
        {
            distance_to_lkgl = calculate_distance(d_last_known_good_lat, d_last_known_good_lon, new_lat, new_lon);

            if (distance_to_lkgl < d_pos_jump_recovery)
                {
                    // Reset jump check when the receiver is back to the last known good location
                    reset_pos_jump_check(new_lat, new_lon, new_alt);
                }
            else
                {
                    d_score.position_jump_score = 2;
                    d_spoofer_score = d_score.total_score();
                    DLOG(INFO) << "POS_JUMP: Spoofer score: " << d_spoofer_score << " - Jump distance: " << jump_distance;

                    // Set last known good location to old coordinates
                    if (d_update_lkgl)
                        {
                            set_last_known_good_location(d_old_lat, d_old_lon, d_old_alt);
                            d_update_lkgl = false;
                        }
                }
        }
    else
        {
            if (d_score.position_jump_score == 2)
                {
                    distance_to_lkgl = calculate_distance(d_last_known_good_lat, d_last_known_good_lon, new_lat, new_lon);
                    if (distance_to_lkgl < d_pos_jump_recovery)
                        {
                            // Reset jump check when the receiver is back to the last known good location
                            reset_pos_jump_check(new_lat, new_lon, new_alt);
                            DLOG(INFO) << "POS_JUMP: Received location within geo-fence, resetting static position score";
                        }
                }
            else
                {
                    set_last_known_good_location(new_lat, new_lon, new_alt);
                }
        }

    set_old_location(new_lat, new_lon, new_alt);
    DLOG(INFO) << "POS_JUMP: Old location updated to: " << new_lat << ", " << new_lon << ", " << new_alt;
    DLOG(INFO) << "POS_JUMP: Last known good location: " << d_last_known_good_lat << ", " << d_last_known_good_lon;
    DLOG(INFO) << "POS_JUMP: Distance to LKGL: " << distance_to_lkgl;
}

// ####### General functions
void SpoofingDetector::set_last_known_good_location(double new_lat, double new_lon, double new_alt)
{
    d_last_known_good_lat = new_lat;
    d_last_known_good_lon = new_lon;
    d_last_known_good_alt = new_alt;

    DLOG(INFO) << "Last known good location updated to: " << new_lat << ", " << new_lon << ", " << new_alt;
}

void SpoofingDetector::set_old_location(double new_lat, double new_lon, double new_alt)
{
    d_old_lat = new_lat;
    d_old_lon = new_lon;
    d_old_alt = new_alt;

    DLOG(INFO) << "Old location updated to: " << new_lat << ", " << new_lon << ", " << new_alt;
}

void SpoofingDetector::reset_pos_jump_check(double new_lat, double new_lon, double new_alt)
{
    d_score.position_jump_score = 0;

    // Set last known good location to current coordinates
    set_last_known_good_location(new_lat, new_lon, new_alt);

    d_score.position_jump_score = 0;
    d_spoofer_score = d_score.total_score();

    d_update_lkgl = true;
}

int SpoofingDetector::get_spoofer_score()
{
    int d_spoofer_score = d_score.total_score();
    DLOG(INFO) << "Total spoofer score: " << d_spoofer_score;
    return d_spoofer_score;
}

long double SpoofingDetector::to_radians(double degree)
{
    // Convert degrees to radians
    long double one_deg = (M_PI) / 180;
    return (one_deg * degree);
}

long double SpoofingDetector::calculate_distance(double lat1, double lon1, double lat2, double lon2)
{
    // Calculate distance between l1 and l2 using Haversine formula
    lat1 = to_radians(lat1);
    lon1 = to_radians(lon1);
    lat2 = to_radians(lat2);
    lon2 = to_radians(lon2);

    // Haversine Formula
    long double dlong = lon2 - lon1;
    long double dlat = lat2 - lat1;

    long double distance = pow(sin(dlat / 2), 2) +
                           cos(lat1) * cos(lat2) *
                               pow(sin(dlong / 2), 2);

    distance = 2 * asin(sqrt(distance));

    long double R = 6371000;

    // Calculate the result
    distance = distance * R;
    return distance;
}