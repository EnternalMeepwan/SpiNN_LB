/*
 * Copyright (c) 2017-2019 The University of Manchester
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//! imports
#include "spin1_api.h"
#include "common-typedefs.h"
#include <data_specification.h>
#include <simulation.h>
#include <debug.h>
#include <circular_buffer.h>
#include <recording.h>

/*
* Data needed by the LBM
*/
#define Q 9

// constant about the speed of lattice
#define cVel 0.57735026919

// velocity vector
int ex[Q] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
int ey[Q] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };

// fi
float fi[Q];
// equilibrium distribution function
float feq[Q];
// post distribution function
float fi_star[Q];

// wi
float wi[Q] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };

// macro fluid density rho init as 1.0
float rho = 1.0;

// u_x u_y, velocity in two dimension
float u_x = 0.0001;
float u_y = 1.0;

// parameter
float nu = 0.000128;
float tau = 0.500384; //(1 + nu * 6) / 2;

/*********************************************************/

/*! multicast routing keys to communicating with neighbours */
uint my_key;

/*! the size of the received buffer at this moment*/
uint current_buffer_size;

/*! buffer used to store spikes */
static circular_buffer input_buffer;
static uint32_t current_key;
static uint32_t current_payload;

// un-matched packets
uint unmatched_packets = 0;

//
uint32_t vertex_index = 0;

// ! position data
int x_pos = 0;
int y_pos = 0;

//! control value, which says how many timer ticks to run for before exiting
static uint32_t simulation_ticks = 0;
static uint32_t time = 0;
data_specification_metadata_t* data = NULL;

//! The recording flags
static uint32_t recording_flags = 0;

//! int as a bool to represent if this simulation should run forever
static uint32_t infinite_run;

// value for turning on and off interrupts
uint cpsr = 0;

//! receive pacekt chekcing
uint received_packets[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

//! human readable definitions of each region in SDRAM
typedef enum regions_e {
    SYSTEM_REGION,
    TRANSMISSIONS,
    POSITION,
    NEIGHBOUR_KEYS,
    VELOCITY,
    VERTEX_INDEX,
    RECORDED_DATA
} regions_e;

//! values for the priority for each callback
typedef enum callback_priorities {
    MC_PACKET = -1,
    SDP = 1,
    TIMER = 2,
    DMA = 3
} callback_priorities;

//! definitions of each element in the transmission region
typedef struct transmission_region {
    uint32_t has_key;
    uint32_t my_key;
} transmission_region_t;

//! definitions of the position of this lattice
typedef struct position_data {
    uint32_t x_position;
    uint32_t y_position;
} position_data_t;

// timer offset
uint32_t timer_offset_wei;

//
typedef struct vertex_index_data {
    uint32_t index;
    uint32_t time_offset;
} vertex_index_data_t;

// diraction = ["me", "n", "w", "s", "e", "nw", "sw", "se", "ne"]
// keys of its neighbours
typedef struct neighbour_keys_region {
    uint32_t n_key;
    uint32_t w_key;
    uint32_t s_key;
    uint32_t e_key;
    uint32_t nw_key;
    uint32_t sw_key;
    uint32_t se_key;
    uint32_t ne_key;
    uint32_t key_mask;
} neighbour_keys_t;
neighbour_keys_t neighbour_keys;

// velocity data region, two dimension
typedef struct velocity_region {
    float U_X;
    float U_Y;
} velocity_t;
velocity_t velocities;

// A union function that convert a int number into float
static inline float int_to_float(int data)
{
    union {
        float x;
        uint32_t y;
    } cast_union;
    cast_union.y = data;
    return cast_union.x;
}

/*
 * Check if the lattice receive the right number of packets and from the right directions
 */
void do_safety_check(void)
{

    cpsr = spin1_int_disable();
    log_info("the unmatched_packets is %d", unmatched_packets);
    if (!(received_packets[0] && received_packets[1] && received_packets[2] && received_packets[3] && received_packets[4] && received_packets[5] && received_packets[6] && received_packets[7])) {
        log_error("didn't receive the correct number of fi_star at time %d, my time offset my timer offset is %d", time, timer_offset_wei);
        rt_error(RTE_SWERR);
    }

    if (unmatched_packets != 56) {
        log_error("did not receive the correct number of un_matched fi_star at time %d", time);
        rt_error(RTE_SWERR);
    }
    spin1_mode_restore(cpsr);
}

/*
 * Receive the fi_star
 */
void read_input_buffer(void)
{
    // wait until receive the correct number of packets
    while (circular_buffer_size(input_buffer) < 128) {
        current_buffer_size = circular_buffer_size(input_buffer);
        for (uint32_t counter = 0; counter < 10000; counter++) {
            //do nothing
        }
        spin1_delay_us(1);
    }
    cpsr = spin1_int_disable();

    for (uint32_t counter = 0; counter < 64; counter++) {
        bool success_get_key = circular_buffer_get_next(input_buffer, &current_key);
        bool success_get_payload = circular_buffer_get_next(input_buffer, &current_payload);
        // get the direction
        uint32_t direction = current_key & 0x00000007;
        // reduce to the base key
        uint32_t base_key = current_key & neighbour_keys.key_mask;
        if (success_get_key && success_get_payload) {
            // diraction = ["me", "n", "w", "s", "e", "nw", "sw", "se", "ne"]
            fi[0] = fi_star[0];
            if (base_key == neighbour_keys.n_key && direction == 0) {
                fi[1] = int_to_float(current_payload);
                received_packets[0] = 1;
            }
            else if (base_key == neighbour_keys.w_key && direction == 1) {
                fi[2] = int_to_float(current_payload);
                received_packets[1] = 1;
            }
            else if (base_key == neighbour_keys.s_key && direction == 2) {
                fi[3] = int_to_float(current_payload);
                received_packets[2] = 1;
            }
            else if (base_key == neighbour_keys.e_key && direction == 3) {
                fi[4] = int_to_float(current_payload);
                received_packets[3] = 1;
            }
            else if (base_key == neighbour_keys.nw_key && direction == 4) {
                fi[5] = int_to_float(current_payload);
                received_packets[4] = 1;
            }
            else if (base_key == neighbour_keys.sw_key && direction == 5) {
                fi[6] = int_to_float(current_payload);
                received_packets[5] = 1;
            }
            else if (base_key == neighbour_keys.se_key && direction == 6) {
                fi[7] = int_to_float(current_payload);
                received_packets[6] = 1;
            }
            else if (base_key == neighbour_keys.ne_key && direction == 7) {
                fi[8] = int_to_float(current_payload);
                received_packets[7] = 1;
            }
            else {
                unmatched_packets += 1;
            }
        }
        else {
            if (!success_get_key)
                log_info("I did not get the key");
            if (!success_get_payload)
                log_info("I did not get the payload");
        }
    }
    spin1_mode_restore(cpsr);
}

void send_with_masked_key()
{
    // diraction = ["n", "w", "s", "e", "nw", "sw", "se", "ne"]
    uint32_t north = float_to_int(fi_star[1]);
    uint32_t west = float_to_int(fi_star[2]);
    uint32_t south = float_to_int(fi_star[3]);
    uint32_t east = float_to_int(fi_star[4]);
    uint32_t northwest = float_to_int(fi_star[5]);
    uint32_t southwest = float_to_int(fi_star[6]);
    uint32_t southeast = float_to_int(fi_star[7]);
    uint32_t northeast = float_to_int(fi_star[8]);

    uint32_t delay = vertex_index % 400;

    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key, north, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 1, west, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 2, south, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 3, east, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 4, northwest, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 5, southwest, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 6, southeast, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 7, northeast, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
}
/*
* Send the fi
*/
void send_state(void)
{
    // reset for next iteration
    for (uint i = 0; i < 8; i++) {
        received_packets[i] = 0;
    }
    unmatched_packets = 0;
    // diraction = ["n", "w", "s", "e", "nw", "sw", "se", "ne"]
    uint delay = vertex_index % 100;
    // add a random delay here
    spin1_delay_us(delay);
    send_with_masked_key();
//    log_info("sent my state via multicast at %d", time);
}

/*
* LBM code start
*/

//! init the distribution function f_i
void initFi(float* fi, float rho)
{
    float A = 1 / (cVel * cVel);
    float B = 1 / 2 * cVel * cVel * cVel * cVel;
    for (int k = 0; k < Q; k++) {
        fi[k] = wi[k] * rho * (1.0 + A * (ex[k] * u_x + ey[k] * u_y) //  A c_ia u_a =  A (c_ix u_x + c_iy u_y)
                                  + B * ((ex[k] * ex[k] - 1.0 / 3.0) // Q_ixx = (c_ix c_ix - 1/3 )
                                            + (ex[k] * ey[k]) // Q_ixy = (c_ix c_iy)
                                            + (ey[k] * ex[k]) // Q_iyx = (c_iy c_ix)
                                            + (ey[k] * ey[k] - 1.0 / 3.0)) //Q_iyy = (c_iy c_iy - 1/3)
                                  );
    }
}

/*
* Compute the equilibrium distribution function f_a^eq for every lattice
* int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
* int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
* float **u_x, float **u_y
*/
void computeFeq(float rho, float u_x, float u_y, float* feq)
{
    float f1 = 3.0;
    float f2 = 9.0 / 2.0;
    float f3 = 3.0 / 2.0;
    float uuDot, euDot, si;

    uuDot = u_x * u_x + u_y * u_y;
    for (int k = 0; k < Q; k++) {
        euDot = ex[k] * u_x + ey[k] * u_y;
        si = f1 * euDot + f2 * euDot * euDot - f3 * uuDot; // cVel^2 = 1/3
        feq[k] = (1.0 + si) * wi[k] * rho;
    }
}

/*
 * compute the macroscopic density rho and velocity u, u = (ux, uy) for every lattice
*/
void computeRho_N_U(float* rho, float* u_x, float* u_y, float* fi)
{
    *rho = 0.0;
    *u_x = 0.0;
    *u_y = 0.0;
    for (int k = 0; k < Q; k++) {
        *rho += fi[k];
        *u_x += ex[k] * fi[k];
        *u_y += ey[k] * fi[k];
    }
    (*u_x) /= (*rho);
    (*u_y) /= (*rho);
}

/*
* ! Collide Step : calculate the update distribution function fi = fi_star - (fi_star - feq) / tau
*/
void collideStep(float* fi, float* feq, float* fi_star, float tau)
{
    for (int k = 0; k < Q; k++) {
        fi_star[k] = fi[k] - (fi[k] - feq[k]) / tau; // fi_star is the post
    }
}

void streamStep()
{
    // stream step start
    send_state();
    read_input_buffer();
    // stream step end
}

float calWholeDensity(float* fi)
{
    float mass = 0.0;
    for (int k = 0; k < Q; k++) {
        mass += fi[k];
    }
    return mass;
}

float calMomentum(float rho, float u_x)
{
    return rho * u_x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

/****f* receive_data
 *
 * SUMMARY
 *  This function is used as a callback for packet received events.
 * receives data from 8 neighbours and updates the fi_star
 *
 * SYNOPSIS
 *  void receive_data (uint key, uint payload)
 *
 * INPUTS
 *   uint key: packet routing key - provided by the RTS
 *   uint payload: packet payload - provided by the RTS
 *
 * SOURCE
 */
void receive_data(uint key, uint payload)
{
    use(key);
    if (!circular_buffer_add(input_buffer, key)) {
        log_error("Could not add key");
    }
    if (!circular_buffer_add(input_buffer, payload)) {
        log_error("Could not add payload");
    }
}

/****f* update
 *
 * SUMMARY
 *
 * SYNOPSIS
 *  void update (uint ticks, uint b)
 *
 * SOURCE
 */
void update(uint ticks, uint b)
{
    use(b);
    use(ticks);

    time++;

//    log_info("on tick %d of %d", time, simulation_ticks);

    // check that the run time hasn't already elapsed and thus needs to be
    // killed
    if ((infinite_run != TRUE) && (time >= simulation_ticks)) {
        // fall into the pause resume mode of operating
        simulation_handle_pause_resume(NULL);

        // Finalise any recordings that are in progress, writing back the final
        // amounts of samples recorded to SDRAM
        if (recording_flags > 0) {
            log_info("updating recording regions");
            recording_record(0, &(u_x), sizeof(float));
            recording_record(0, &(u_y), sizeof(float));
            recording_do_timestep_update(time);
            recording_finalise();
        }

        log_info("Simulation complete.");

        // switch to state where host is ready to read
        simulation_ready_to_read();

        return;
    }
    //    float momentum, mass;

    if (time == 0) {
        initFi(fi, rho);
        computeFeq(rho, u_x, u_y, feq);
    }
    else {
        do_safety_check();
    }
    computeRho_N_U(&rho, &u_x, &u_y, fi);
    computeFeq(rho, u_x, u_y, feq);
    collideStep(fi, feq, fi_star, tau);
    streamStep();
}



//! \brief this method is to catch strange behaviour
//! \param[in] key: the key being received
//! \param[in] unknown: second arg with no state. set to zero by default
void receive_data_void(uint key, uint unknown)
{
    use(key);
    use(unknown);
    log_error("this should never ever be done");
}

static bool initialize(uint32_t* timer_period, uint32_t* timer_offset)
{
    log_info("Initialise: started");

    // Get the address this core's DTCM data starts at from SRAM
    data = data_specification_get_data_address();

    // Read the header
    if (!data_specification_read_header(data)) {
        log_error("failed to read the data spec header");
        return false;
    }

    // Get the timing details and set up the simulation interface
    if (!simulation_initialise(
            data_specification_get_region(SYSTEM_REGION, data),
            APPLICATION_NAME_HASH, timer_period, &simulation_ticks,
            &infinite_run, &time, SDP, DMA)) {
        return false;
    }

    // initialise transmission keys
    transmission_region_t* transmission_sdram = data_specification_get_region(TRANSMISSIONS, data);
    if (!transmission_sdram->has_key) {
        log_error(
            "this lattice cell can't affect anything, deduced as an error,"
            "please fix the application fabric and try again");
        return false;
    }
    my_key = transmission_sdram->my_key;
    log_info("my key is %d", my_key);

    // read test data for init tick
    position_data_t* position_data_sdram = data_specification_get_region(POSITION, data);
    x_pos = position_data_sdram->x_position;
    y_pos = position_data_sdram->y_position;
    log_info("I'm in position ( %d, %d)", x_pos, y_pos);

    neighbour_keys_t* neighbour_data_sdram = data_specification_get_region(NEIGHBOUR_KEYS, data);
    neighbour_keys = *neighbour_data_sdram;

    log_info("my n key is %d", neighbour_keys.n_key);

    velocity_t* velocity_sdram = data_specification_get_region(VELOCITY, data);
    u_x = velocity_sdram->U_X;
    u_y = velocity_sdram->U_Y;
    log_info("my U_X = %f", u_x);
    log_info("my U_Y = %f", u_y);

    vertex_index_data_t* vertex_index_sdram = data_specification_get_region(VERTEX_INDEX, data);
    vertex_index = vertex_index_sdram->index;
    *timer_offset = vertex_index_sdram->time_offset;
    log_info("vertex_index = %d", vertex_index);



    // initialise my input_buffer for receiving packets
    input_buffer = circular_buffer_initialize(2048);
    if (input_buffer == 0) {
        return false;
    }
    log_info("input_buffer initialised");

    void* recording_region = data_specification_get_region(RECORDED_DATA, data);
    bool success = recording_initialize(recording_region, &recording_flags);
    log_info("Recording flags = 0x%08x", recording_flags);
    return success;
}

/****f* c_main
 *
 * SUMMARY
 *  This function is called at application start-up.
 *  It is used to register event callbacks and begin the simulation.
 *
 * SYNOPSIS
 *  int c_main()
 *
 * SOURCE
 */
void c_main(void)
{
    log_info("starting lattice cell");

    // Load DTCM data
    uint32_t timer_period, timer_offset;

    // initialise the model
    if (!initialize(&timer_period, &timer_offset)) {
        log_error("Error in initialisation - exiting!");
        rt_error(RTE_SWERR);
    }

    // set timer tick value to configured value
    log_info("setting timer to execute every %d microseconds", timer_period);
    log_info("setting timer offset as %d microsends", timer_offset);
    spin1_set_timer_tick_and_phase(timer_period, timer_offset);

    // register callbacks
    spin1_callback_on(MCPL_PACKET_RECEIVED, receive_data, MC_PACKET);
    spin1_callback_on(TIMER_TICK, update, TIMER);
    timer_offset_wei = timer_offset;
    // start execution
    log_info("Starting\n");

    // Start the time at "-1" so that the first tick will be 0
    time = UINT32_MAX;

    simulation_run();
}/*
 * Copyright (c) 2017-2019 The University of Manchester
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//! imports
#include "spin1_api.h"
#include "common-typedefs.h"
#include <data_specification.h>
#include <simulation.h>
#include <debug.h>
#include <circular_buffer.h>
#include <recording.h>

/*
* Data needed by the LBM
*/
#define Q 9

// constant about the speed of lattice
#define cVel 0.57735026919

// velocity vector
int ex[Q] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
int ey[Q] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };

// fi
float fi[Q];
// equilibrium distribution function
float feq[Q];
// post distribution function
float fi_star[Q];

// wi
float wi[Q] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };

// macro fluid density rho init as 1.0
float rho = 1.0;

// u_x u_y, velocity in two dimension
float u_x = 0.0001;
float u_y = 1.0;

// parameter
float nu = 0.000128;
float tau = 0.500384; //(1 + nu * 6) / 2;

/*********************************************************/

/*! multicast routing keys to communicating with neighbours */
uint my_key;

/*! the size of the received buffer at this moment*/
uint current_buffer_size;

/*! buffer used to store spikes */
static circular_buffer input_buffer;
static uint32_t current_key;
static uint32_t current_payload;

// un-matched packets
uint unmatched_packets = 0;

//
uint32_t vertex_index = 0;

// ! position data
int x_pos = 0;
int y_pos = 0;

//! control value, which says how many timer ticks to run for before exiting
static uint32_t simulation_ticks = 0;
static uint32_t time = 0;
data_specification_metadata_t* data = NULL;

//! The recording flags
static uint32_t recording_flags = 0;

//! int as a bool to represent if this simulation should run forever
static uint32_t infinite_run;

// value for turning on and off interrupts
uint cpsr = 0;

//! receive pacekt chekcing
uint received_packets[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

//! human readable definitions of each region in SDRAM
typedef enum regions_e {
    SYSTEM_REGION,
    TRANSMISSIONS,
    POSITION,
    NEIGHBOUR_KEYS,
    VELOCITY,
    VERTEX_INDEX,
    RECORDED_DATA
} regions_e;

//! values for the priority for each callback
typedef enum callback_priorities {
    MC_PACKET = -1,
    SDP = 1,
    TIMER = 2,
    DMA = 3
} callback_priorities;

//! definitions of each element in the transmission region
typedef struct transmission_region {
    uint32_t has_key;
    uint32_t my_key;
} transmission_region_t;

//! definitions of the position of this lattice
typedef struct position_data {
    uint32_t x_position;
    uint32_t y_position;
} position_data_t;

// timer offset
uint32_t timer_offset_wei;

//
typedef struct vertex_index_data {
    uint32_t index;
    uint32_t time_offset;
} vertex_index_data_t;

// diraction = ["me", "n", "w", "s", "e", "nw", "sw", "se", "ne"]
// keys of its neighbours
typedef struct neighbour_keys_region {
    uint32_t n_key;
    uint32_t w_key;
    uint32_t s_key;
    uint32_t e_key;
    uint32_t nw_key;
    uint32_t sw_key;
    uint32_t se_key;
    uint32_t ne_key;
    uint32_t key_mask;
} neighbour_keys_t;
neighbour_keys_t neighbour_keys;

// velocity data region, two dimension
typedef struct velocity_region {
    float U_X;
    float U_Y;
} velocity_t;
velocity_t velocities;

// A union function that convert a int number into float
static inline float int_to_float(int data)
{
    union {
        float x;
        uint32_t y;
    } cast_union;
    cast_union.y = data;
    return cast_union.x;
}

/*
 * Check if the lattice receive the right number of packets and from the right directions
 */
void do_safety_check(void)
{

    cpsr = spin1_int_disable();
    log_info("the unmatched_packets is %d", unmatched_packets);
    if (!(received_packets[0] && received_packets[1] && received_packets[2] && received_packets[3] && received_packets[4] && received_packets[5] && received_packets[6] && received_packets[7])) {
        log_error("didn't receive the correct number of fi_star at time %d, my time offset my timer offset is %d", time, timer_offset_wei);
        rt_error(RTE_SWERR);
    }

    if (unmatched_packets != 56) {
        log_error("did not receive the correct number of un_matched fi_star at time %d", time);
        rt_error(RTE_SWERR);
    }
    spin1_mode_restore(cpsr);
}

/*
 * Receive the fi_star
 */
void read_input_buffer(void)
{
    // wait until receive the correct number of packets
    while (circular_buffer_size(input_buffer) < 128) {
        current_buffer_size = circular_buffer_size(input_buffer);
        for (uint32_t counter = 0; counter < 10000; counter++) {
            //do nothing
        }
        spin1_delay_us(1);
    }
    cpsr = spin1_int_disable();

    for (uint32_t counter = 0; counter < 64; counter++) {
        bool success_get_key = circular_buffer_get_next(input_buffer, &current_key);
        bool success_get_payload = circular_buffer_get_next(input_buffer, &current_payload);
        // get the direction
        uint32_t direction = current_key & 0x00000007;
        // reduce to the base key
        uint32_t base_key = current_key & neighbour_keys.key_mask;
        if (success_get_key && success_get_payload) {
            // diraction = ["me", "n", "w", "s", "e", "nw", "sw", "se", "ne"]
            fi[0] = fi_star[0];
            if (base_key == neighbour_keys.n_key && direction == 0) {
                fi[1] = int_to_float(current_payload);
                received_packets[0] = 1;
            }
            else if (base_key == neighbour_keys.w_key && direction == 1) {
                fi[2] = int_to_float(current_payload);
                received_packets[1] = 1;
            }
            else if (base_key == neighbour_keys.s_key && direction == 2) {
                fi[3] = int_to_float(current_payload);
                received_packets[2] = 1;
            }
            else if (base_key == neighbour_keys.e_key && direction == 3) {
                fi[4] = int_to_float(current_payload);
                received_packets[3] = 1;
            }
            else if (base_key == neighbour_keys.nw_key && direction == 4) {
                fi[5] = int_to_float(current_payload);
                received_packets[4] = 1;
            }
            else if (base_key == neighbour_keys.sw_key && direction == 5) {
                fi[6] = int_to_float(current_payload);
                received_packets[5] = 1;
            }
            else if (base_key == neighbour_keys.se_key && direction == 6) {
                fi[7] = int_to_float(current_payload);
                received_packets[6] = 1;
            }
            else if (base_key == neighbour_keys.ne_key && direction == 7) {
                fi[8] = int_to_float(current_payload);
                received_packets[7] = 1;
            }
            else {
                unmatched_packets += 1;
            }
        }
        else {
            if (!success_get_key)
                log_info("I did not get the key");
            if (!success_get_payload)
                log_info("I did not get the payload");
        }
    }
    spin1_mode_restore(cpsr);
}

void send_with_masked_key()
{
    // diraction = ["n", "w", "s", "e", "nw", "sw", "se", "ne"]
    uint32_t north = float_to_int(fi_star[1]);
    uint32_t west = float_to_int(fi_star[2]);
    uint32_t south = float_to_int(fi_star[3]);
    uint32_t east = float_to_int(fi_star[4]);
    uint32_t northwest = float_to_int(fi_star[5]);
    uint32_t southwest = float_to_int(fi_star[6]);
    uint32_t southeast = float_to_int(fi_star[7]);
    uint32_t northeast = float_to_int(fi_star[8]);

    uint32_t delay = vertex_index % 400;

    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key, north, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 1, west, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 2, south, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 3, east, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 4, northwest, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 5, southwest, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 6, southeast, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
    spin1_delay_us(delay);
    while (!spin1_send_mc_packet(my_key + 7, northeast, WITH_PAYLOAD)) {
        spin1_delay_us(1);
    }
}
/*
* Send the fi
*/
void send_state(void)
{
    // reset for next iteration
    for (uint i = 0; i < 8; i++) {
        received_packets[i] = 0;
    }
    unmatched_packets = 0;
    // diraction = ["n", "w", "s", "e", "nw", "sw", "se", "ne"]
    uint delay = vertex_index % 100;
    // add a random delay here
    spin1_delay_us(delay);
    send_with_masked_key();
//    log_info("sent my state via multicast at %d", time);
}

/*
* LBM code start
*/

//! init the distribution function f_i
void initFi(float* fi, float rho)
{
    float A = 1 / (cVel * cVel);
    float B = 1 / 2 * cVel * cVel * cVel * cVel;
    for (int k = 0; k < Q; k++) {
        fi[k] = wi[k] * rho * (1.0 + A * (ex[k] * u_x + ey[k] * u_y) //  A c_ia u_a =  A (c_ix u_x + c_iy u_y)
                                  + B * ((ex[k] * ex[k] - 1.0 / 3.0) // Q_ixx = (c_ix c_ix - 1/3 )
                                            + (ex[k] * ey[k]) // Q_ixy = (c_ix c_iy)
                                            + (ey[k] * ex[k]) // Q_iyx = (c_iy c_ix)
                                            + (ey[k] * ey[k] - 1.0 / 3.0)) //Q_iyy = (c_iy c_iy - 1/3)
                                  );
    }
}

/*
* Compute the equilibrium distribution function f_a^eq for every lattice
* int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
* int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
* float **u_x, float **u_y
*/
void computeFeq(float rho, float u_x, float u_y, float* feq)
{
    float f1 = 3.0;
    float f2 = 9.0 / 2.0;
    float f3 = 3.0 / 2.0;
    float uuDot, euDot, si;

    uuDot = u_x * u_x + u_y * u_y;
    for (int k = 0; k < Q; k++) {
        euDot = ex[k] * u_x + ey[k] * u_y;
        si = f1 * euDot + f2 * euDot * euDot - f3 * uuDot; // cVel^2 = 1/3
        feq[k] = (1.0 + si) * wi[k] * rho;
    }
}

/*
 * compute the macroscopic density rho and velocity u, u = (ux, uy) for every lattice
*/
void computeRho_N_U(float* rho, float* u_x, float* u_y, float* fi)
{
    *rho = 0.0;
    *u_x = 0.0;
    *u_y = 0.0;
    for (int k = 0; k < Q; k++) {
        *rho += fi[k];
        *u_x += ex[k] * fi[k];
        *u_y += ey[k] * fi[k];
    }
    (*u_x) /= (*rho);
    (*u_y) /= (*rho);
}

/*
* ! Collide Step : calculate the update distribution function fi = fi_star - (fi_star - feq) / tau
*/
void collideStep(float* fi, float* feq, float* fi_star, float tau)
{
    for (int k = 0; k < Q; k++) {
        fi_star[k] = fi[k] - (fi[k] - feq[k]) / tau; // fi_star is the post
    }
}

void streamStep()
{
    // stream step start
    send_state();
    read_input_buffer();
    // stream step end
}

float calWholeDensity(float* fi)
{
    float mass = 0.0;
    for (int k = 0; k < Q; k++) {
        mass += fi[k];
    }
    return mass;
}

float calMomentum(float rho, float u_x)
{
    return rho * u_x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

/****f* receive_data
 *
 * SUMMARY
 *  This function is used as a callback for packet received events.
 * receives data from 8 neighbours and updates the fi_star
 *
 * SYNOPSIS
 *  void receive_data (uint key, uint payload)
 *
 * INPUTS
 *   uint key: packet routing key - provided by the RTS
 *   uint payload: packet payload - provided by the RTS
 *
 * SOURCE
 */
void receive_data(uint key, uint payload)
{
    use(key);
    if (!circular_buffer_add(input_buffer, key)) {
        log_error("Could not add key");
    }
    if (!circular_buffer_add(input_buffer, payload)) {
        log_error("Could not add payload");
    }
}

/****f* update
 *
 * SUMMARY
 *
 * SYNOPSIS
 *  void update (uint ticks, uint b)
 *
 * SOURCE
 */
void update(uint ticks, uint b)
{
    use(b);
    use(ticks);

    time++;

//    log_info("on tick %d of %d", time, simulation_ticks);

    // check that the run time hasn't already elapsed and thus needs to be
    // killed
    if ((infinite_run != TRUE) && (time >= simulation_ticks)) {
        // fall into the pause resume mode of operating
        simulation_handle_pause_resume(NULL);

        // Finalise any recordings that are in progress, writing back the final
        // amounts of samples recorded to SDRAM
        if (recording_flags > 0) {
            log_info("updating recording regions");
            recording_record(0, &(u_x), sizeof(float));
            recording_record(0, &(u_y), sizeof(float));
            recording_do_timestep_update(time);
            recording_finalise();
        }

        log_info("Simulation complete.");

        // switch to state where host is ready to read
        simulation_ready_to_read();

        return;
    }
    //    float momentum, mass;

    if (time == 0) {
        initFi(fi, rho);
        computeFeq(rho, u_x, u_y, feq);
    }
    else {
        do_safety_check();
    }
    computeRho_N_U(&rho, &u_x, &u_y, fi);
    computeFeq(rho, u_x, u_y, feq);
    collideStep(fi, feq, fi_star, tau);
    streamStep();
}



//! \brief this method is to catch strange behaviour
//! \param[in] key: the key being received
//! \param[in] unknown: second arg with no state. set to zero by default
void receive_data_void(uint key, uint unknown)
{
    use(key);
    use(unknown);
    log_error("this should never ever be done");
}

static bool initialize(uint32_t* timer_period, uint32_t* timer_offset)
{
    log_info("Initialise: started");

    // Get the address this core's DTCM data starts at from SRAM
    data = data_specification_get_data_address();

    // Read the header
    if (!data_specification_read_header(data)) {
        log_error("failed to read the data spec header");
        return false;
    }

    // Get the timing details and set up the simulation interface
    if (!simulation_initialise(
            data_specification_get_region(SYSTEM_REGION, data),
            APPLICATION_NAME_HASH, timer_period, &simulation_ticks,
            &infinite_run, &time, SDP, DMA)) {
        return false;
    }

    // initialise transmission keys
    transmission_region_t* transmission_sdram = data_specification_get_region(TRANSMISSIONS, data);
    if (!transmission_sdram->has_key) {
        log_error(
            "this lattice cell can't affect anything, deduced as an error,"
            "please fix the application fabric and try again");
        return false;
    }
    my_key = transmission_sdram->my_key;
    log_info("my key is %d", my_key);

    // read test data for init tick
    position_data_t* position_data_sdram = data_specification_get_region(POSITION, data);
    x_pos = position_data_sdram->x_position;
    y_pos = position_data_sdram->y_position;
    log_info("I'm in position ( %d, %d)", x_pos, y_pos);

    neighbour_keys_t* neighbour_data_sdram = data_specification_get_region(NEIGHBOUR_KEYS, data);
    neighbour_keys = *neighbour_data_sdram;

    log_info("my n key is %d", neighbour_keys.n_key);

    velocity_t* velocity_sdram = data_specification_get_region(VELOCITY, data);
    u_x = velocity_sdram->U_X;
    u_y = velocity_sdram->U_Y;
    log_info("my U_X = %f", u_x);
    log_info("my U_Y = %f", u_y);

    vertex_index_data_t* vertex_index_sdram = data_specification_get_region(VERTEX_INDEX, data);
    vertex_index = vertex_index_sdram->index;
    *timer_offset = vertex_index_sdram->time_offset;
    log_info("vertex_index = %d", vertex_index);



    // initialise my input_buffer for receiving packets
    input_buffer = circular_buffer_initialize(2048);
    if (input_buffer == 0) {
        return false;
    }
    log_info("input_buffer initialised");

    void* recording_region = data_specification_get_region(RECORDED_DATA, data);
    bool success = recording_initialize(recording_region, &recording_flags);
    log_info("Recording flags = 0x%08x", recording_flags);
    return success;
}

/****f* c_main
 *
 * SUMMARY
 *  This function is called at application start-up.
 *  It is used to register event callbacks and begin the simulation.
 *
 * SYNOPSIS
 *  int c_main()
 *
 * SOURCE
 */
void c_main(void)
{
    log_info("starting lattice cell");

    // Load DTCM data
    uint32_t timer_period, timer_offset;

    // initialise the model
    if (!initialize(&timer_period, &timer_offset)) {
        log_error("Error in initialisation - exiting!");
        rt_error(RTE_SWERR);
    }

    // set timer tick value to configured value
    log_info("setting timer to execute every %d microseconds", timer_period);
    log_info("setting timer offset as %d microsends", timer_offset);
    spin1_set_timer_tick_and_phase(timer_period, timer_offset);

    // register callbacks
    spin1_callback_on(MCPL_PACKET_RECEIVED, receive_data, MC_PACKET);
    spin1_callback_on(TIMER_TICK, update, TIMER);
    timer_offset_wei = timer_offset;
    // start execution
    log_info("Starting\n");

    // Start the time at "-1" so that the first tick will be 0
    time = UINT32_MAX;

    simulation_run();
}
