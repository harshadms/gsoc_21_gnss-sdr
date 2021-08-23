#include "spoofing_detector.h"
#include "channel_event.h"

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

#include <cmath>  // for floor, fmod, rint, ceil
#include <map>
#include <numeric>

SpoofingDetectorConfig::SpoofingDetectorConfig()
{
    // ####### Position consistency check variables

    max_jump_distance = 100;   // meters
    geo_fence_radius = 15;     // meters
    velocity_difference = 15;  // m/s
    pos_error_threshold = 10;  // meters
    min_altitude = -10;        // meters
    max_altitude = 20000;      // meters
    min_ground_speed = 0;      // m/s
    max_ground_speed = 200;    // m/s
    static_lat = 0;            //degrees
    static_lon = 0;            //degrees
    static_alt = 0;            //meters

    position_check = false;
    dump_pvt_checks_results = false;
    static_pos_check = false;
    use_aux_peak = false;
    enable_apt = false;

    clk_offset_vector_size = 1000;
    clk_offset_error = 20;

    std::string filename = "./pos_consistency_results.mat";
}

SpoofingDetector::SpoofingDetector()
{
}

SpoofingDetector::SpoofingDetector(const SpoofingDetectorConfig* conf_)
{
    DLOG(INFO) << "Spoofing detector for PVT initialized";

    d_max_jump_distance = conf_->max_jump_distance;
    d_geo_fence_radius = conf_->geo_fence_radius;
    d_velocity_difference = conf_->velocity_difference;
    d_pos_error_threshold = conf_->pos_error_threshold;

    // Set predetermined location
    d_static_pvt.lat = conf_->static_lat;
    d_static_pvt.lon = conf_->static_lon;
    d_static_pvt.alt = conf_->static_alt;

    d_dump_pvt_checks_results = conf_->dump_pvt_checks_results;
    d_position_check = conf_->position_check;
    d_static_pos_check = conf_->static_pos_check;

    // Position abnormalities check paramaters
    d_min_altitude = conf_->min_altitude;
    d_max_altitude = conf_->max_altitude;
    d_min_ground_speed = conf_->min_ground_speed;
    d_max_ground_speed = conf_->max_ground_speed;

    d_clk_offset_vector_size = conf_->clk_offset_vector_size;
    d_clk_offset_error = conf_->clk_offset_error;

    d_cno_threshold = conf_->cno_threshold;

    d_spoofer_score = 0;
    d_checked_velocity_pairs = 0;  // Number of velocity measurements checked. Velocity measurements are processed in pairs (old and new)

    // Used to set old coordinates and last known good location
    d_first_record = true;

    // Used to decide whether to update LKGL
    d_update_lkgl = true;

    boost::posix_time::ptime d_pvt_epoch(boost::gregorian::date(1970, 1, 1));
}

// ####### Position consistency functions
void SpoofingDetector::check_PVT_consistency()
{
    // Public function
    // PVT block calls this function. Actual consistency checks are triggered from here
    if (d_static_pos_check)
        {
            if (validate_location_proximity(&d_new_pvt, &d_static_pvt, 1))
                {
                    if (d_score.static_pos_check_score > 0)
                        {
                            DLOG(INFO) << "STATIC_POS: Received location within geo-fence, resetting static position score";
                        }

                    d_score.static_pos_check_score = 0;
                }
        }

    check_position_jump();

    check_velocity_consistency();

    abnormal_position_checks();

    check_time();

    if (d_first_record)
        {
            d_first_record = false;
        }
    update_old_pvt();
    d_spoofer_score = d_score.total_score();
}

void SpoofingDetector::abnormal_position_checks()
{
    // A collection of model based position abnormalities check Min/Max altitude and speed
    d_score.abnormal_position_score = 0;

    //if (d_new_pvt.alt < d_min_altitude) d_score.abnormal_position_score += 0.25;
    //if (d_new_pvt.alt > d_max_altitude) d_score.abnormal_position_score += 0.25;
    if (d_new_pvt.speed_over_ground < d_min_ground_speed) d_score.abnormal_position_score += 0.25;
    if (d_new_pvt.speed_over_ground > d_max_ground_speed) d_score.abnormal_position_score += 0.25;

    DLOG(INFO) << "ABNORMAL_CHECK: " << d_score.abnormal_position_score;
}

void SpoofingDetector::check_velocity_consistency()
{
    /// Compares the reported velocity with the position pairs. If the projected coordinates do not match the received coordinate velocity error is increased
    ++d_checked_velocity_pairs;

    if (!check_propagated_pos())
        {
            ++d_velocity_error;
        }
    d_score.velocity_check_score = d_velocity_error / d_checked_velocity_pairs;
    DLOG(INFO) << "VELOCITY_COMPARE: Confidence " << d_velocity_error << "/" << d_checked_velocity_pairs;
}

bool SpoofingDetector::validate_location_proximity(const PvtSol* pvtsol1, const PvtSol* pvtsol2, int test_id)
{
    long double distance = calculate_distance(pvtsol1, pvtsol2);

    switch (test_id)
        {
        case 1:  // Validate location's proximity to the configured static location
            {
                DLOG(INFO) << "STATIC_POS: static_coords (" << pvtsol1->lat << "," << pvtsol1->lon << ") received (" << pvtsol2->lat << "," << pvtsol2->lon << ") Distance: " << distance;

                if (distance < d_geo_fence_radius)
                    {
                        d_score.static_pos_check_score = 0;
                        return true;
                    }
                else
                    {
                        d_score.static_pos_check_score = distance;
                        return false;
                    }
            }
        case 2:  // Compare the distance between the old and the new solutions with configured max jump distance
            {
                DLOG(INFO) << "POS_JUMP," << pvtsol1->lat << "," << pvtsol1->lon << ") received ("
                           << pvtsol2->lat << "," << pvtsol2->lon << ") Distance: "
                           << distance << " Spoofer score: " << d_score.position_jump_score;

                if (distance < d_max_jump_distance)
                    {
                        d_score.position_jump_score = 0;
                        return true;
                    }
                else
                    {
                        d_score.position_jump_score = distance;
                        return false;
                    }
            }

        case 3:  // Check if the calculated solution is in the proximity of the lkg solution
            {
                return distance < d_pos_error_threshold;
            }
        }

    return false;
}

bool SpoofingDetector::check_position_jump()
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
                    DLOG(INFO) << "POS_JUMP: location update: " << d_new_pvt.lat << ", " << d_new_pvt.lon << ", " << d_new_pvt.alt;
                }

            d_first_record = false;
            return false;
        }

    if (d_score.position_jump_score > 0)
        {
            if (validate_location_proximity(&d_lkg_pvt, &d_new_pvt, 3))
                {
                    // Reset jump check when the receiver is back to the last known good location
                    reset_pos_jump_check();
                }
        }
    else
        {
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
                                    DLOG(INFO) << "POS_JUMP: Spoofer score: " << d_score.position_jump_score;

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
            update_lkg_pvt(false);
        }

    return false;
}

// Propagates the PVT solution using the reported velocity and validates the location proximity
bool SpoofingDetector::check_propagated_pos()
{
    PvtSol temp_pvt;
    double dt = (d_new_pvt.tstamp - d_old_pvt.tstamp) / 1000;

    double metersPerDegLat = 111111.0;
    double metersPerRadLat = metersPerDegLat * 180 / M_PI;
    double metersPerDegLon = metersPerDegLat * cos(d_old_pvt.lat);
    double metersPerRadLon = metersPerDegLon * 180 / M_PI;

    temp_pvt.lat = d_old_pvt.lat + d_old_pvt.vel_y * dt / metersPerRadLat;
    temp_pvt.lon = d_old_pvt.lon + d_old_pvt.vel_x * dt / metersPerRadLon;
    temp_pvt.alt = d_old_pvt.alt + d_old_pvt.vel_z * dt;

    double distance_error = calculate_distance(&temp_pvt, &d_new_pvt);

    DLOG(INFO) << "PROPAGATE_POS: Pro " << temp_pvt.lat << "," << temp_pvt.lon << ", " << temp_pvt.alt;
    DLOG(INFO) << "PROPAGATE_POS: Recv " << d_new_pvt.lat << "," << d_new_pvt.lon << ", " << d_new_pvt.alt;
    DLOG(INFO) << "PROPAGATE_POS: Error: " << distance_error;

    return validate_location_proximity(&temp_pvt, &d_new_pvt, 3);
}

// Check calculated UTC time
void SpoofingDetector::check_time()
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
bool SpoofingDetector::check_clock_offset(double clk_offset, double clk_drift)
{
    ClockOffset offset;
    offset.offset = clk_offset * 1e9;
    offset.drift = clk_drift * 1e-3;  //aging per nanosec

    offset.timestamp = SpoofingDetector::CurrentTime_nanoseconds();

    d_clock_offsets_vector.push_back(offset);

    if (d_clock_offsets_vector.size() < d_clk_offset_vector_size)
        {
            d_clock_offsets_vector.push_back(offset);
            return false;
        }

    if (d_clock_offsets_vector.size() > d_clk_offset_vector_size)
        {
            int no_elements_to_remove = d_clock_offsets_vector.size() - d_clk_offset_vector_size;
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
    bool spoofing_flag = false;

    // Implement a moving average of offset error
    d_clock_offest_errors_vector.push_back(offsetError);
    if (d_clock_offest_errors_vector.size() > d_clk_offset_vector_size / 10)
        {
            int no_elements_to_remove = d_clock_offest_errors_vector.size() - d_clk_offset_vector_size / 10;
            d_clock_offest_errors_vector.erase(d_clock_offest_errors_vector.begin(), d_clock_offest_errors_vector.begin() + no_elements_to_remove);

            double avg_error = std::accumulate(d_clock_offest_errors_vector.begin(), d_clock_offest_errors_vector.end(), 0) / d_clock_offest_errors_vector.size();

            if (avg_error > d_clk_offset_error)
                {
                    spoofing_flag = true;
                    d_score.clk_offset = spoofing_flag;
                    DLOG(INFO) << "CLK_OFFSET: Spoofing detected at " << SpoofingDetector::CurrentTime_nanoseconds() << " ns  Recv offset: " << clk_offset * 1e9 << " - projected: " << offset_propd << " - error: " << offsetError;
                }
        }
    DLOG(INFO) << "CLK_PLT," << SpoofingDetector::CurrentTime_nanoseconds() << "," << clk_offset * 1e9 << "," << offset_propd << "," << offsetError << "," << spoofing_flag;
    return spoofing_flag;
}

// Clock jump
double SpoofingDetector::check_RX_clock()
{
    if (d_first_record)
        {
            DLOG(INFO) << "TLM_CHECKS: Enabled clock jump check";
            set_old_clock();
            set_lkg_clock(false);
            d_first_record = false;
            return false;
        }

    return check_clock_jump();
}


double SpoofingDetector::check_clock_jump()
{
    double sample_diff = d_new_clock.sample_counter - d_lkg_clock.sample_counter;
    double sample_diff_time = (sample_diff * 1000) / d_gnss_synchro->fs;

    uint32_t tow_diff = d_new_clock.tow - d_lkg_clock.tow;

    uint32_t time_diff = abs(tow_diff - sample_diff_time);
    bool flag = false;

    if (time_diff > 20)
        {
            flag = true;
            d_gnss_synchro->Clock_jump = time_diff;
            set_old_clock();
            set_lkg_clock(true);
            return time_diff;
        }
    DLOG(INFO) << "TLM_CHECKS_" << d_gnss_synchro->PRN << "," << sample_diff_time << "," << tow_diff << "," << flag;
    set_old_clock();
    set_lkg_clock(false);
    return 0;
}


void SpoofingDetector::update_clock_info(uint64_t sample_counter, uint32_t tow, uint32_t wn)
{
    d_new_clock.sample_counter = sample_counter;
    d_new_clock.tow = tow;
    d_new_clock.wn = wn;
    //DLOG(INFO) << "TLM_CHECKS: Received new clock info - Channel " << channel_id << " : PRN " << PRN << " = " << sample_counter << " " << tow;
}


void SpoofingDetector::set_old_clock()
{
    d_old_clock.sample_counter = d_new_clock.sample_counter;
    d_old_clock.tow = d_new_clock.tow;
    d_old_clock.wn = d_new_clock.wn;
}


void SpoofingDetector::set_lkg_clock(bool set_old)
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

// CNO Check
void SpoofingDetector::check_CNO(std::vector<double> cno_vector)
{
    double sum = 0.0;
    double mean = 0.0;
    double standardDeviation = 0.0;

    int i;

    for (i = 0; i < cno_vector.size(); ++i)
        {
            sum += cno_vector[i];
        }

    mean = sum / 10;

    for (i = 0; i < cno_vector.size(); ++i)
        standardDeviation += pow(cno_vector[i] - mean, 2);

    standardDeviation = sqrt(standardDeviation / 10);

    if (standardDeviation < d_cno_threshold)
        {
            DLOG(INFO) << "CNO_STD: STD Exceeds threshold - " << standardDeviation;
            d_score.cno = standardDeviation;
        }
}


// ####### General functions
void SpoofingDetector::update_pvt(double lat, double lon, double alt,
    double vel_x, double vel_y, double vel_z,
    double speed_over_ground, double heading,
    uint32_t tstamp, boost::posix_time::ptime utc_time)
{
    d_new_pvt.lat = lat;
    d_new_pvt.lon = lon;
    d_new_pvt.alt = alt;

    d_new_pvt.vel_x = vel_x;
    d_new_pvt.vel_y = vel_y;
    d_new_pvt.vel_z = vel_z;

    d_new_pvt.speed_over_ground = speed_over_ground;
    d_new_pvt.heading = heading;

    d_new_pvt.tstamp = tstamp;
    d_new_pvt.utc_time = utc_time;

    DLOG(INFO)
        << "New PVT: " << lat << ", " << lon << ", " << alt
        << ", " << vel_x << ", " << vel_y << ", " << vel_z
        << ", " << speed_over_ground << ", " << heading << ", " << tstamp;
}

void SpoofingDetector::update_old_pvt()
{
    d_old_pvt.lat = d_new_pvt.lat;
    d_old_pvt.lon = d_new_pvt.lon;
    d_old_pvt.alt = d_new_pvt.alt;
    d_old_pvt.vel_x = d_new_pvt.vel_x;
    d_old_pvt.vel_y = d_new_pvt.vel_y;
    d_old_pvt.vel_z = d_new_pvt.vel_z;
    d_old_pvt.speed_over_ground = d_new_pvt.speed_over_ground;
    d_old_pvt.heading = d_new_pvt.heading;
    d_old_pvt.tstamp = d_new_pvt.tstamp;
    d_old_pvt.utc_time = d_new_pvt.utc_time;

    DLOG(INFO) << "Old pvt updated to: " << d_old_pvt.lat << ", " << d_old_pvt.lon << ", " << d_old_pvt.alt;
}

void SpoofingDetector::update_lkg_pvt(bool set_old)
{
    if (set_old)
        {
            d_lkg_pvt.lat = d_old_pvt.lat;
            d_lkg_pvt.lon = d_old_pvt.lon;
            d_lkg_pvt.alt = d_old_pvt.alt;
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
            d_lkg_pvt.lat = d_new_pvt.lat;
            d_lkg_pvt.lon = d_new_pvt.lon;
            d_lkg_pvt.alt = d_new_pvt.alt;
            d_lkg_pvt.vel_x = d_new_pvt.vel_x;
            d_lkg_pvt.vel_y = d_new_pvt.vel_y;
            d_lkg_pvt.vel_z = d_new_pvt.vel_z;
            d_lkg_pvt.speed_over_ground = d_new_pvt.speed_over_ground;
            d_lkg_pvt.heading = d_new_pvt.heading;
            d_lkg_pvt.tstamp = d_new_pvt.tstamp;
            d_lkg_pvt.utc_time = d_new_pvt.utc_time;
        }


    DLOG(INFO) << "LKG updated to: " << d_lkg_pvt.lat << ", " << d_lkg_pvt.lon << ", " << d_lkg_pvt.alt;
}

void SpoofingDetector::reset_pos_jump_check()
{
    // Set last known good location to current coordinates
    update_lkg_pvt(false);  // Set old = false, hence set new location as lkg

    d_score.position_jump_score = 0;
    d_spoofer_score = d_score.total_score();
    d_update_lkgl = true;

    DLOG(INFO) << "POS_JUMP: Received location within geo-fence, resetting static position score";
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

long double SpoofingDetector::calculate_distance_ECEF(const PvtSol* pvtsol1, const PvtSol* pvtsol2)
{
    return sqrt(pow((pvtsol1->lon - pvtsol2->lon), 2) + pow((pvtsol1->lat - pvtsol2->lat), 2) + pow((pvtsol1->alt - pvtsol2->alt), 2));
}

long double SpoofingDetector::calculate_distance(const PvtSol* pvtsol1, const PvtSol* pvtsol2)
{
    // Calculate distance between l1 and l2 using Haversine formula
    long double lat1 = to_radians(pvtsol1->lat);
    long double lon1 = to_radians(pvtsol1->lon);

    long double lat2 = to_radians(pvtsol2->lat);
    long double lon2 = to_radians(pvtsol2->lon);

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

void SpoofingDetector::set_msg_queue(std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> control_queue)
{
    DLOG(INFO) << "Setting message queue for spoofing detector object";
    control_queue_ = std::move(control_queue);
}

void SpoofingDetector::stop_tracking(int channel)
{
    control_queue_->push(pmt::make_any(channel_event_make(channel, 3)));
}

//#shift everything here, use synchro tracking counter and TOW at current symbol for clock jump