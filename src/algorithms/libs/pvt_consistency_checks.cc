/*!
 * \file pvt_consistency_checks.cc
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

PVTConsistencyChecks::PVTConsistencyChecks()
{
}

PVTConsistencyChecks::PVTConsistencyChecks(const PVTConsistencyChecksConf* conf_)
{
    DLOG(INFO) << "Spoofing detector for PVT initialized";

    d_max_jump_distance = conf_->max_jump_distance;
    d_geo_fence_radius = conf_->geo_fence_radius;
    d_velocity_difference = conf_->velocity_difference;
    d_pos_error_threshold = conf_->pos_error_threshold;

    // Set predetermined location
    d_static_pvt.ecef_y = conf_->static_lat;
    d_static_pvt.ecef_x = conf_->static_lon;
    d_static_pvt.ecef_z = conf_->static_alt;

    d_dump_pvt_checks_results = conf_->dump_pvt_checks_results;
    d_position_check = conf_->position_check;
    d_static_pos_check = conf_->static_pos_check;

    // Position abnormalities check paramaters
    d_min_altitude = conf_->min_altitude;
    d_max_altitude = conf_->max_altitude;
    d_min_ground_speed = conf_->min_ground_speed;
    d_max_ground_speed = conf_->max_ground_speed;

    d_spoofer_score = 0;
    d_checked_velocity_pairs = 0;  // Number of velocity measurements checked. Velocity measurements are processed in pairs (old and new)

    // Used to set old coordinates and last known good location
    d_first_record = true;

    // Used to decide whether to update LKGL
    d_update_lkgl = true;


    boost::posix_time::ptime d_pvt_epoch(boost::gregorian::date(1970, 1, 1));
}

// ####### Position consistency functions
void PVTConsistencyChecks::check_PVT_consistency()
{
    // Public function
    // PVT block calls this function. Actual consistency checks are triggered from here
    if (d_static_pos_check)
        {
            if (validate_location_proximity(&d_new_pvt, &d_static_pvt, 1))
                {
                    if (d_score.static_pos_check_score == 1)
                        {
                            DLOG(INFO) << "STATIC_POS: Received location within geo-fence, resetting static position score";
                        }

                    d_score.static_pos_check_score = 0;
                }
            else
                {
                    d_score.static_pos_check_score = 1;
                    d_spoofer_score = d_score.total_score();
                }
        }

    if (check_position_jump())
        {
            d_score.position_jump_score = 2;
        }

    check_velocity_consistency();

    abnormal_position_checks();

    check_time();

    if (d_first_record)
        {
            d_first_record = false;
        }

    d_spoofer_score = d_score.total_score();
}

void PVTConsistencyChecks::abnormal_position_checks()
{
    // A collection of model based position abnormalities check Min/Max altitude and speed
    d_score.abnormal_position_score = 0;

    //if (d_new_pvt.ecef_z < d_min_altitude) d_score.abnormal_position_score += 0.25;
    //if (d_new_pvt.ecef_z > d_max_altitude) d_score.abnormal_position_score += 0.25;
    if (d_new_pvt.speed_over_ground < d_min_ground_speed) d_score.abnormal_position_score += 0.25;
    if (d_new_pvt.speed_over_ground > d_max_ground_speed) d_score.abnormal_position_score += 0.25;

    DLOG(INFO) << "ABNORMAL_CHECK: " << d_score.abnormal_position_score;
}

void PVTConsistencyChecks::check_velocity_consistency()
{
    /// Compares the reported velocity with the position pairs. If the projected coordinates do not match the received coordinate velocity error is increased
    ++d_checked_velocity_pairs;

    if (!check_propagated_pos())
        {
            ++d_velocity_error;
        }
    d_score.velocity_check_score = d_velocity_error / d_checked_velocity_pairs;
    DLOG(INFO) << "VELOCITY_COMPARE: Confidence " << d_velocity_error << "/" << d_checked_velocity_pairs;
    update_old_pvt();
}

bool PVTConsistencyChecks::validate_location_proximity(const PvtSol* pvtsol1, const PvtSol* pvtsol2, int test_id)
{
    long double distance = calculate_distance_ECEF(pvtsol1, pvtsol2);

    switch (test_id)
        {
        case 1:  // Validate location's proximity to the configured static location
            {
                DLOG(INFO) << "STATIC_POS: static_coords (" << pvtsol1->ecef_y << "," << pvtsol1->ecef_x << ") received (" << pvtsol2->ecef_y << "," << pvtsol2->ecef_x << ") Distance: " << distance;
                return distance < d_geo_fence_radius;
            }
        case 2:  // Compare the distance between the old and the new solutions with configured max jump distance
            {
                DLOG(INFO) << "POS_JUMP: Old (" << pvtsol1->ecef_y << "," << pvtsol1->ecef_x << ") received ("
                           << pvtsol2->ecef_y << "," << pvtsol2->ecef_x << ") Distance: "
                           << distance << " Spoofer score: " << get_spoofer_score();

                return distance < d_max_jump_distance;
            }

        case 3:  // Check if the calculated solution is in the proximity of the lkg solution
            {
                return distance < d_pos_error_threshold;
            }
        }
}

bool PVTConsistencyChecks::check_position_jump()
{
    bool is_spoofing = false;

    if (d_first_record)
        {
            update_old_pvt();

            // In case of a cold start there is no way to know if the received location is legit. Hence a naive way is to check for spoofer score.
            // The score will be 0 if no other technique detects spoofing
            if (get_spoofer_score() == 0)
                {
                    // Set last known good location to current coordinates
                    update_lkg_pvt(false);  // Set old = false, hence set new location as lkg
                    DLOG(INFO) << "POS_JUMP: location update: " << d_new_pvt.ecef_y << ", " << d_new_pvt.ecef_x << ", " << d_new_pvt.ecef_z;
                }

            d_first_record = false;
            return false;
        }


    if (!validate_location_proximity(&d_old_pvt, &d_new_pvt, 2))
        {
            if (validate_location_proximity(&d_lkg_pvt, &d_new_pvt, 3))
                {
                    // Reset jump check when the receiver is back to the last known good location
                    reset_pos_jump_check();
                }
            else
                {
                    // A naive way of checking if the jump is legitimate jump or is caused by spoofing.
                    // (A legitimate jump will occur if the receiver looses lock and re-acquires it)
                    //uint32_t dt = d_new_pvt.tstamp - d_lkg_pvt.tstamp;

                    if (!check_propagated_pos())
                        {
                            is_spoofing = true;
                        }

                    if (is_spoofing)
                        {
                            DLOG(INFO) << "POS_JUMP: Spoofer score: " << d_spoofer_score;

                            // Set last known good location to old coordinates
                            if (d_update_lkgl)
                                {
                                    update_lkg_pvt(true);  // Set old = true, hence set old location as lkg
                                    d_update_lkgl = false;
                                }
                            return true;
                        }
                }
        }
    else
        {
            if (d_score.position_jump_score == 2)
                {
                    if (validate_location_proximity(&d_lkg_pvt, &d_new_pvt, 3))
                        {
                            // Reset jump check when the receiver is back to the last known good location
                            reset_pos_jump_check();
                            DLOG(INFO) << "POS_JUMP: Received location within geo-fence, resetting static position score";
                        }
                }
            else
                {
                    update_lkg_pvt(false);
                }
        }

    return false;
}

// Propagates the PVT solution using the reported velocity and validates the location proximity
bool PVTConsistencyChecks::check_propagated_pos()
{
    PvtSol temp_pvt;
    double dt = (d_new_pvt.tstamp - d_old_pvt.tstamp) / 1000;

    temp_pvt.ecef_y = d_old_pvt.ecef_y + d_old_pvt.vel_x * dt;  // / metersPerRadLat;
    temp_pvt.ecef_x = d_old_pvt.ecef_x + d_old_pvt.vel_y * dt;  // / metersPerRadLon;
    temp_pvt.ecef_z = d_old_pvt.ecef_z + d_old_pvt.vel_z * dt;

    double distance_error = calculate_distance_ECEF(&temp_pvt, &d_new_pvt);

    DLOG(INFO) << "PROPAGATE_POS: Pro " << temp_pvt.ecef_y << "," << temp_pvt.ecef_x << ", " << temp_pvt.ecef_z;
    DLOG(INFO) << "PROPAGATE_POS: Recv " << d_new_pvt.ecef_y << "," << d_new_pvt.ecef_x << ", " << d_new_pvt.ecef_z;
    DLOG(INFO) << "PROPAGATE_POS: Error: " << distance_error;

    return validate_location_proximity(&temp_pvt, &d_new_pvt, 3);
}

void PVTConsistencyChecks::check_time()
{
    boost::posix_time::ptime now(boost::posix_time::microsec_clock::universal_time());

    if ((now - d_new_pvt.utc_time).total_seconds() < -18)
        {
            DLOG(INFO) << "UTC_TIME_CHECK: Calculated UTC time is " << (now - d_new_pvt.utc_time).total_seconds() << " in future";
        }
    else if ((now - d_new_pvt.utc_time).total_seconds() > 18)
        {
            DLOG(INFO) << "UTC_TIME_CHECK: Calculated UTC time is " << (now - d_new_pvt.utc_time).total_seconds() << " in past";
        }
}

// Clock offset
void PVTConsistencyChecks::check_clock_offset(double clk_offset, double clk_drift)
{
    ClockOffset offset;
    offset.offset = clk_offset * 1e9;
    offset.drift = clk_drift * 1e-3;  //aging per sec

    offset.timestamp = PVTConsistencyChecks::CurrentTime_nanoseconds();

    d_clock_offsets_vector.push_back(offset);

    // TO:DO add a config param for this value
    if (d_clock_offsets_vector.size() < 2500)
        {
            d_clock_offsets_vector.push_back(offset);
            return;
        }

    if (d_clock_offsets_vector.size() > 2500)
        {
            int no_elements_to_remove = d_clock_offsets_vector.size() - 2500;
            d_clock_offsets_vector.erase(d_clock_offsets_vector.begin(), d_clock_offsets_vector.begin() + no_elements_to_remove);
        }

    // Credits - PNT-Integrity library
    double driftExp = 0.0;
    double driftVar = 0.0;
    unsigned int i = 0;

    for (auto it = d_clock_offsets_vector.begin(); it != (d_clock_offsets_vector.end() - 1); ++it)
        {
            driftExp += it->drift;
            driftVar += pow(it->drift, 2);
            ++i;
        }

    driftExp = driftExp / i;
    driftVar = driftVar / i;
    driftVar = driftVar - pow(driftExp, 2);
    /// \note driftVar - pow(driftExp,2) can sometimes be slightly negative
    /// due to quantization, set all negative values to 0
    if (driftVar < 0)
        {
            driftVar = 0.0;
        }

    double dt = ((d_clock_offsets_vector.rbegin())->timestamp - (d_clock_offsets_vector.rbegin() + 1)->timestamp);
    double offset_propd = (d_clock_offsets_vector.rbegin() + 1)->offset + (driftExp * dt);
    double offsetError = fabs(offset_propd - d_clock_offsets_vector.rbegin()->offset);

    DLOG(INFO) << "CLK_OFFSET: Recv offset: " << clk_offset * 1e9 << " - projected: " << offset_propd << " - error: " << offsetError;
}

// ####### General functions
void PVTConsistencyChecks::update_pvt(const std::array<double, 3>& pos,
    const std::array<double, 3>& vel,
    double speed_over_ground, double heading,
    uint32_t tstamp, boost::posix_time::ptime utc_time)
{
    d_new_pvt.ecef_y = pos[1];
    d_new_pvt.ecef_x = pos[0];
    d_new_pvt.ecef_z = pos[2];

    d_new_pvt.vel_x = vel[0];
    d_new_pvt.vel_y = vel[1];
    d_new_pvt.vel_z = vel[2];

    d_new_pvt.speed_over_ground = sqrt(pow(vel[0], 2) + pow(vel[1], 2) + pow(vel[1], 2));
    d_new_pvt.heading = 180 + 180 / M_PI * (atan2(-vel[0], -vel[2]));
    d_new_pvt.tstamp = tstamp;
    d_new_pvt.utc_time = utc_time;

    DLOG(INFO)
        << "New PVT: " << pos[1] << ", " << pos[0] << ", " << pos[2]
        << ", " << vel[0] << ", " << vel[1] << ", " << vel[2]
        << ", " << speed_over_ground << ", " << heading << ", " << tstamp;
}

void PVTConsistencyChecks::update_old_pvt()
{
    d_old_pvt.ecef_y = d_new_pvt.ecef_y;
    d_old_pvt.ecef_x = d_new_pvt.ecef_x;
    d_old_pvt.ecef_z = d_new_pvt.ecef_z;
    d_old_pvt.vel_x = d_new_pvt.vel_x;
    d_old_pvt.vel_y = d_new_pvt.vel_y;
    d_old_pvt.vel_z = d_new_pvt.vel_z;
    d_old_pvt.speed_over_ground = d_new_pvt.speed_over_ground;
    d_old_pvt.heading = d_new_pvt.heading;
    d_old_pvt.tstamp = d_new_pvt.tstamp;
    d_old_pvt.utc_time = d_new_pvt.utc_time;

    DLOG(INFO) << "Old pvt updated to: " << d_old_pvt.ecef_y << ", " << d_old_pvt.ecef_x << ", " << d_old_pvt.ecef_z;
}

void PVTConsistencyChecks::update_lkg_pvt(bool set_old)
{
    if (set_old)
        {
            d_lkg_pvt.ecef_y = d_old_pvt.ecef_y;
            d_lkg_pvt.ecef_x = d_old_pvt.ecef_x;
            d_lkg_pvt.ecef_z = d_old_pvt.ecef_z;
            d_lkg_pvt.vel_x = d_old_pvt.vel_x;
            d_lkg_pvt.vel_y = d_old_pvt.vel_y;
            d_lkg_pvt.vel_z = d_old_pvt.vel_z;
            d_lkg_pvt.speed_over_ground = d_old_pvt.speed_over_ground;
            d_lkg_pvt.heading = d_old_pvt.heading;
            d_lkg_pvt.tstamp = d_old_pvt.tstamp;
            d_lkg_pvt.utc_time = d_old_pvt.utc_time;
        }
    else
        {
            d_lkg_pvt.ecef_y = d_new_pvt.ecef_y;
            d_lkg_pvt.ecef_x = d_new_pvt.ecef_x;
            d_lkg_pvt.ecef_z = d_new_pvt.ecef_z;
            d_lkg_pvt.vel_x = d_new_pvt.vel_x;
            d_lkg_pvt.vel_y = d_new_pvt.vel_y;
            d_lkg_pvt.vel_z = d_new_pvt.vel_z;
            d_lkg_pvt.speed_over_ground = d_new_pvt.speed_over_ground;
            d_lkg_pvt.heading = d_new_pvt.heading;
            d_lkg_pvt.tstamp = d_new_pvt.tstamp;
            d_lkg_pvt.utc_time = d_new_pvt.utc_time;
        }


    DLOG(INFO) << "LKG updated to: " << d_lkg_pvt.ecef_y << ", " << d_lkg_pvt.ecef_x << ", " << d_lkg_pvt.ecef_z;
}

void PVTConsistencyChecks::reset_pos_jump_check()
{
    // Set last known good location to current coordinates
    update_lkg_pvt(false);  // Set old = false, hence set new location as lkg

    d_score.position_jump_score = 0;
    d_spoofer_score = d_score.total_score();

    d_update_lkgl = true;
}

int PVTConsistencyChecks::get_spoofer_score()
{
    int d_spoofer_score = d_score.total_score();
    DLOG(INFO) << "Total spoofer score: " << d_spoofer_score;
    return d_spoofer_score;
}

long double PVTConsistencyChecks::to_radians(double degree)
{
    // Convert degrees to radians
    long double one_deg = (M_PI) / 180;
    return (one_deg * degree);
}

long double PVTConsistencyChecks::calculate_distance_ECEF(const PvtSol* pvtsol1, const PvtSol* pvtsol2)
{
    return sqrt(pow((pvtsol1->ecef_x - pvtsol2->ecef_x), 2) + pow((pvtsol1->ecef_y - pvtsol2->ecef_y), 2) + pow((pvtsol1->ecef_z - pvtsol2->ecef_z), 2));
}