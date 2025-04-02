#ifndef MAGNET2D
#define MAGNET2D

#include "my_vectors.h"
#include "geometry.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
  #define M_PI 3.141592653589793238462643383279502984
#endif

typedef struct magnet_2d {

  struct polygon core;
  struct polygon wing;
  double B_core;
  double B_wing;

} magnet_2d;

typedef struct configuration {

  int num_magnets;
  struct magnet_2d *magnets;

} configuration;


struct magnet_2d make_magnet_2d(int num_core, const double *core_x, const double *core_y, const double *core_z,
                              int num_wing, const double *wing_x, const double *wing_y, const double *wing_z,
                              double B_goal, int type);

struct configuration make_empty_configuration(void);

void add_magnet_to_configuration(struct magnet_2d *mag, struct configuration *config);

int is_inside_magnet(const struct state_vector *point, const struct magnet_2d *mag);

// struct vector compute_B_from_configuration(const struct vector *point, const struct configuration *config);

struct configuration make_configuration(int num_magnets, const int *num_core, const int *num_wing,
                                        const double *x_core, const double *y_core, const double *z_core,
                                        const double *x_wing, const double *y_wing, const double *z_wing,
                                        const double *B_goal, const int *types);

#endif