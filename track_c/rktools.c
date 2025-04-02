#include "rktools.h"

void track_particles(double *ret_data,
                     int num_particles,
                     double h,
                     double T_max,
                     int subsampling_rate,
                     const struct track_limits *limits,
                     const struct configuration *config,
                     const struct particle *prtcls,
                     int *num_steps){

    // auxiliary variables
    int i, j;

    // the number of samples calculated
    int num_t = (int) (T_max/h);

    // the number of stored instances
    int num_t_stored = (int) (num_t/subsampling_rate);

    // printf("T_max = %e\n", T_max);
    // printf("h = %e\n", h);
    // printf("subsampling_rate = %d\n", subsampling_rate);
    // printf("Number of times to store = %d\n", num_t_stored);

    // make the return data
    // double *ret_data = (double*) calloc(6*num_particles*num_t_stored, sizeof(double));

    #pragma omp parallel
    {
        // make the track object
        struct track tr = make_empty_track();

        // allocate the memory for the tracks
        free(tr.states);
        tr.num_steps = num_t;
        tr.states = (struct state_vector*) calloc(num_t, sizeof(struct state_vector));

        #pragma omp for
		for(i = 0; i < num_particles; ++i){

            track_particle(&tr, 0.0, h, T_max, limits, config, &prtcls[i]);

            #pragma omp critical
            {
                for(j = 0; j < num_t_stored; ++j){
                    
                    ret_data[6*num_t_stored*i + 6*j] = tr.states[subsampling_rate*j].x;
                    ret_data[6*num_t_stored*i + 6*j + 1] = tr.states[subsampling_rate*j].y;
                    ret_data[6*num_t_stored*i + 6*j + 2] = tr.states[subsampling_rate*j].z;
                    ret_data[6*num_t_stored*i + 6*j + 3] = tr.states[subsampling_rate*j].px;
                    ret_data[6*num_t_stored*i + 6*j + 4] = tr.states[subsampling_rate*j].py;
                    ret_data[6*num_t_stored*i + 6*j + 5] = tr.states[subsampling_rate*j].pz;
                    

                }

                num_steps[i] = (int) (tr.num_steps/subsampling_rate);
                
            }
        }

        free(tr.states);
        tr.states = NULL;
    }
    
    
    return;
}

void track_particle(struct track *tr,
                    double t_0,
                    double h,
                    double T_max,
                    const struct track_limits *limits,
                    const struct configuration *config,
                    const struct particle *prtcl){

    // declare all local variables
    int num_steps_max, check_limits, n;
    struct state_vector k_1, k_2, k_3, k_4;
    double t = t_0;

    // determine the number of steps (maximum)
    num_steps_max = (int) (T_max/h);
    
   
    // the first state is the initial one
    tr->states[0] = prtcl->y0;
    
    for (n = 1; n < num_steps_max; ++n){

        // compute the auxiliary variables
        // compute the auxiliary variables
        k_1 = eval_f(t, tr->states[n-1], config, prtcl);
        k_2 = eval_f(t + 0.5*h, add_scaled_state_vectors(1.0, &tr->states[n-1], 0.5*h, &k_1), config, prtcl);
        k_3 = eval_f(t + 0.5*h, add_scaled_state_vectors(1.0, &tr->states[n-1], 0.5*h, &k_2), config, prtcl);
        k_4 = eval_f(t + h, add_scaled_state_vectors(1.0, &tr->states[n-1], h, &k_3), config, prtcl);

        // sum up
        tr->states[n].x = tr->states[n-1].x + h/6.0*(k_1.x + 2.0*k_2.x + 2.0*k_3.x + k_4.x);
        tr->states[n].y = tr->states[n-1].y + h/6.0*(k_1.y + 2.0*k_2.y + 2.0*k_3.y + k_4.y);
        tr->states[n].z = tr->states[n-1].z + h/6.0*(k_1.z + 2.0*k_2.z + 2.0*k_3.z + k_4.z);
        tr->states[n].px = tr->states[n-1].px + h/6.0*(k_1.px + 2.0*k_2.px + 2.0*k_3.px + k_4.px);
        tr->states[n].py = tr->states[n-1].py + h/6.0*(k_1.py + 2.0*k_2.py + 2.0*k_3.py + k_4.py);
        tr->states[n].pz = tr->states[n-1].pz + h/6.0*(k_1.pz + 2.0*k_2.pz + 2.0*k_3.pz + k_4.pz);

        // tr->states[n] = add_scaled_state_vectors(1.0, &tr->states[n-1], h/6.0, &k_1);
        // tr->states[n] = add_scaled_state_vectors(1.0, &tr->states[n], h/3.0, &k_2);
        // tr->states[n] = add_scaled_state_vectors(1.0, &tr->states[n], h/3.0, &k_3);
        // tr->states[n] = add_scaled_state_vectors(1.0, &tr->states[n], h/6.0, &k_4);

        // increment the time step
        t += h;

        // check limits
        check_limits = limits_reached(&tr->states[n], limits);

        // break loop if limits are reached
        if (check_limits == 1){

           break;
        }
    }

    // set the number of steps
    tr->num_steps = n;

    return;

}

int limits_reached(const struct state_vector *Y, const struct track_limits *lim){

    if (Y->x >= lim->x_max){
        return 1;
    }
    if (Y->y >= lim->y_max){
        return 1;
    }
    if (Y->z >= lim->z_max){
        return 1;
    }
    if (Y->x <= lim->x_min){
        return 1;
    }
    if (Y->y <= lim->y_min){
        return 1;
    }
    if (Y->z <= lim->z_min){
        return 1;
    }
    return 0;
}


struct state_vector eval_f(double t, struct state_vector Y, const struct configuration *config, const struct particle *prtcl){

    // local variables
    double c_sq, c_sqsq, m_sq, p_norm_sq, den, fac;
    double Bx = 0.0;
    double By = 0.0;
    double Bz = 0.0;
    int check;
    struct state_vector dY_ret;

    // some auxiliary variables
    c_sq = CONST_C*CONST_C;
    c_sqsq = c_sq*c_sq;
    m_sq = prtcl->m*prtcl->m;
    p_norm_sq = Y.px*Y.px + Y.py*Y.py + Y.pz*Y.pz;
    den = sqrt(m_sq*c_sqsq + p_norm_sq*c_sq);
    fac = c_sq/den;

    // make the evaluation point
    // point = make_vector_3d(Y.x, Y.y, Y.z);

    // compute the B field
    for(int i = 0; i < config->num_magnets; ++i){
        check = is_inside_magnet(&Y, &config->magnets[i]);

        if(check == 1){
            By = config->magnets[i].B_core;
            break;
        }
        if(check == 2){
            By = config->magnets[i].B_wing;
            break;
        }
    }

    // fill the return vector
    dY_ret.x = Y.px*fac;
    dY_ret.y = Y.py*fac;
    dY_ret.z = Y.pz*fac;
    dY_ret.px = prtcl->q*(Y.py*Bz*fac - Y.pz*By*fac);
    dY_ret.py = prtcl->q*(Y.pz*Bx*fac - Y.px*Bz*fac);
    dY_ret.pz = prtcl->q*(Y.px*By*fac - Y.py*Bx*fac);


    return dY_ret;
}

struct particle *make_particles(int num_particles, const double *m, const double *q, const double *y0){

    // allocate memory for the particles
    struct particle *prtcls = (struct particle*) calloc(num_particles, sizeof(struct particle));

    // fill them
    for (int i = 0; i < num_particles; ++i){
        prtcls[i].m = m[i];
        prtcls[i].q = q[i];
        
        prtcls[i].y0.x = y0[6*i];
        prtcls[i].y0.y = y0[6*i+1];
        prtcls[i].y0.z = y0[6*i+2];

        prtcls[i].y0.px = y0[6*i+3];
        prtcls[i].y0.py = y0[6*i+4];
        prtcls[i].y0.pz = y0[6*i+5];

    }

    return prtcls;
}

struct track_limits make_track_limits(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max){

    // make the return data
    struct track_limits lim;

    lim.x_max = x_max;
    lim.x_min = x_min;
    lim.y_max = y_max;
    lim.y_min = y_min;
    lim.z_min = z_min;
    lim.z_max = z_max;

    return lim;

}