#ifndef GEOMETRY
#define GEOMETRY

#include "my_vectors.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
  #define M_PI 3.141592653589793238462643383279502984
#endif

typedef struct polygon {

  int num_corners;
  struct vector *cp;

} polygon;


struct polygon make_empty_polygon(void);

struct polygon make_polygon(int num_corners, const double *cp_x, const double *cp_y, const double *cp_z);

void add_corner_point(struct polygon *poly, double x, double y, double z);

void print_polygon(const struct polygon *poly);

double compute_polygon_surface(const struct polygon *poly);

int is_inside_polygon(const struct state_vector *state, const struct polygon *poly);

#endif