/*
 * This is a basic implementation of a lattice Boltzmann method on SpiNNaker.
 * This project is currently being released under the GPL 3.0 license. Use it Free as in
 * Freedom. 
 *
 * The up-to-date information of SpiNNaker Project can be found here:
 * https://spinnakermanchester.github.io/
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
* Data regions for lattice Boltzmann
*/

//! nine neighbours
#define Q 9

//! constant about the speed of lattice 1/sqrt(3)
#define cVel 0.57735026919

//! velocity vector in two dimensions
int ex[Q] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
int ey[Q] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };

//! fi distribution function
float fi[Q];

//! equilibrium distribution function
float feq[Q];

//! post distribution function
float fi_star[Q];

//! wi weighting factors for nine direction in D2Q9 model
float wi[Q] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };

//! macro fluid density rho init as 1.0
float rho = 1.0;

//! u_x u_y, velocity in two dimension
float u_x = 0.0001;
float u_y = 1.0;

//! parameter you need to change this according to the simulation
float nu = 0.000128; // kinematic viscosity
float tau = 0.500384; // relaxation time (reverse time)

/*
* Data regions for SpiNNaker data specification and communication
*/

/*! multicast routing keys to communicating with neighbours */
uint my_key;

/*! the size of the received buffer at this moment*/
uint current_buffer_size;

/*! buffer used to store spikes */
static circular_buffer input_buffer;
static uint32_t current_key;
static uint32_t current_payload;

//! un-matched packets for each time step
uint unmatched_packets = 0;

//! the vertex index of this lattice
uint32_t vertex_index = 0;

//! position of this lattice
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

//! receive packets checking
uint received_packets[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

//! human readable definitions of each region in SDRAM
typedef enum regions_e {
    SYSTEM_REGION, // system configuration
    TRANSMISSIONS, // communication specification
    POSITION, // (x,y)
    NEIGHBOUR_KEYS, // the routing keys of this eight neighbours
    VELOCITY, // the initialized velocities of this lattice
    VERTEX_INDEX, // the index of this lattice in the SpiNNaker's network
    RECORDED_DATA // recording the result
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

//! timer offset. We add a timer offset to reduce the peak packets rate
uint32_t __timer_offset;


typedef struct vertex_index_data {
    uint32_t index;
    uint32_t time_offset;
} vertex_index_data_t;

//! direction = ["me", "n", "w", "s", "e", "nw", "sw", "se", "ne"]
//! keys of its neighbours
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

//! initialized velocity data region, two dimension
typedef struct velocity_region {
    float U_X;
    float U_Y;
} velocity_t;
velocity_t velocities;

//! A union function that convert a int number into float, the other function that transfer a float into int has been
// implemented in the imports
static inline float int_to_float(int data)
{
    union {
        float x;
        uint32_t y;
    } cast_union;
    cast_union.y = data;
    return cast_union.x;
}

 /****f* do_safety_check
 *
 * SUMMARY
 *  Check if the lattice receive the right number of packets and from the right directions. For every time step
 *  the lattice should receive usefull packet from each of 9 direction and 56 unmatched packets from the 9 directions.
 * 
 * SYNOPSIS
 *  void do_safety_check(void)
 *
 * SOURCE
 */
void do_safety_check(void)
{

    cpsr = spin1_int_disable();

    //! check the usefull packets
    if (!(received_packets[0] && received_packets[1] && received_packets[2] && received_packets[3] && received_packets[4] && received_packets[5] && received_packets[6] && received_packets[7])) {
        log_error("didn't receive the correct number of fi_star at time %d, my time offset my timer offset is %d", time, __timer_offset);
        rt_error(RTE_SWERR);
    }

    if (unmatched_packets != 56) {
    
    //! check the unmatched packets
        log_error("did not receive the correct number of un_matched fi_star at time %d", time);
        rt_error(RTE_SWERR);
    }
    spin1_mode_restore(cpsr);
}


 /****f* read_input_buffer
 *
 * SUMMARY
 * read the post-distribution function from the packets payloads. This function firstly will first wait until
 * there is 128 buffered data in the buffer, where 64 of them are keys and the other 64 of them are payloads.
 * After we get the right number of buffer, it will read one key from the buffer, extrate the direction, if it
 * need that packet, it will extrate the floating point number from the next buffer (the corresponding payload)
 * as the post-distribution function. if the key and direction do not match with its need, the unmatched packets
 * number plus 1 for the safety check later.
 * 
 * 
 * SYNOPSIS
 *  read_input_buffer(void)
 *
 * SOURCE
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

 /****f* send_with_masked_key
 *
 * SUMMARY
 * This function will send the its distribution function via multicast. it will first transfer the floating
 * point number into integer. Then between each send, we introduce a random delay to reduce the peak packets
 * rate.
 * 
 * SYNOPSIS
 *  send_with_masked_key(void)
 *
 * SOURCE
 */
void send_with_masked_key()
{
    // transfer the floating-point distribution function to integer format for sending the packets
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
* Send the fi after a reset the receiving information for last timstep and a random delay
*/
void send_state(void)
{
    // reset for next iteration
    for (uint i = 0; i < 8; i++) {
        received_packets[i] = 0;
    }
    unmatched_packets = 0;

    // add a random delay here
    uint delay = vertex_index % 100;
    spin1_delay_us(delay);

    send_with_masked_key();
    log_info("sent my state via multicast at %d", time);
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
* Compute/initialize the equilibrium distribution function f_a^eq for this lattice
*/
void computeFeq(float rho, float u_x, float u_y, float* feq)
{
    float f1 = 3.0; // first factor
    float f2 = 9.0 / 2.0; // second factor
    float f3 = 3.0 / 2.0; // third factor
    float uuDot, euDot, si;

    uuDot = u_x * u_x + u_y * u_y;
    for (int k = 0; k < Q; k++) {
        euDot = ex[k] * u_x + ey[k] * u_y;
        si = f1 * euDot + f2 * euDot * euDot - f3 * uuDot; // cVel^2 = 1/3
        feq[k] = (1.0 + si) * wi[k] * rho;
    }
}

/*
 * compute the macroscopic density rho and velocity u, u = (ux, uy) for this lattice
 * rho = \sum_0^9 f_i
 * u_x = \sum_0^9 e_x * f_i / rho
 * u_x = \sum_0^9 e_y * f_i / rho
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
        fi_star[k] = fi[k] - (fi[k] - feq[k]) / tau; // fi_star is the post-distribution function
    }
}

//
void streamStep()
{
    // send the distribution function via multicast
    send_state();
    // receive the distribution function of its neighbours as post-distribution function
    read_input_buffer();
}

// calculate the density of this lattice. This is usefull if you want to check the physical correctness
float calWholeDensity(float* fi)
{
    float mass = 0.0;
    for (int k = 0; k < Q; k++) {
        mass += fi[k];
    }
    return mass;
}

// calculate the total momentum of this lattice. This is usefull if to want to check physical correctness
float calMomentum(float rho, float u_x)
{
    return rho * u_x;
}


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
    // add the key of the incoming packet to the buffer
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

   log_info("on tick %d of %d", time, simulation_ticks);

    // check that the run time hasn't already elapsed and thus needs to be
    // killed
    if ((infinite_run != TRUE) && (time >= simulation_ticks)) {
        // fall into the pause resume mode of operating
        simulation_handle_pause_resume(NULL);

        // Finalise any recordings that are in progress, writing back the final
        // amounts of samples recorded to SDRAM
        if (recording_flags > 0) {
            log_info("updating recording regions");
            // We only record the result here for the sake of speed. Though we can record result from every 
            // time step but it would take a big amount of time get the data back to host.
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

    if (time == 0) {
        // we only need to initialize the distribution function, the velocity initialization done at front end
        initFi(fi, rho);
        computeFeq(rho, u_x, u_y, feq);
    }
    else {
        // da a safety check for the last timestep if it do not pass this, The result will def wrong.
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

    // read neighbours' routing keys
    neighbour_keys_t* neighbour_data_sdram = data_specification_get_region(NEIGHBOUR_KEYS, data);
    neighbour_keys = *neighbour_data_sdram;
    log_info("my n key is %d", neighbour_keys.n_key);

    // read the initialized velocity of this lattice
    velocity_t* velocity_sdram = data_specification_get_region(VELOCITY, data);
    u_x = velocity_sdram->U_X;
    u_y = velocity_sdram->U_Y;
    log_info("my U_X = %f", u_x);
    log_info("my U_Y = %f", u_y);

    // read the index of this lattice in the communication network
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
    __timer_offset = timer_offset;
    // start execution
    log_info("Starting\n");

    // Start the time at "-1" so that the first tick will be 0
    time = UINT32_MAX;

    simulation_run();
}
