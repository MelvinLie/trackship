#ifndef RK_TOOLS
#define RK_TOOLS

#include "my_vectors.h"
#include "magnet2D.h"

#include <stdio.h>
#include <omp.h>
#include <math.h>

#define CONST_C 299792458.0

typedef struct track_limits {

  double x_max;
  double x_min;
  double y_min;
  double y_max;
  double z_min;
  double z_max;

} track_limits;

typedef struct particle{

  double m;
  double q;
  struct state_vector y0;

} particle;

void track_particle(struct track *tr,
                    double t_0,
                    double step_size,
                    double T_max,
                    const struct track_limits *limits,
                    const struct configuration *config,
                    const struct particle *prtcl);


int limits_reached(const struct state_vector *Y, const struct track_limits *lim);

struct state_vector eval_f(double t, struct state_vector Y, const struct configuration *config, const struct particle *prtcl);

struct particle *make_particles(int num_particles, const double *m, const double *q, const double *y0);

struct track_limits make_track_limits(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max);

void track_particles(double *ret_data, int num_particles, double h, double T_max, int subsampling_rate, const struct track_limits *limits, const struct configuration *config, const struct particle *prtcls, int *num_samples);

#endif