#ifndef MY_VECTORS
#define MY_VECTORS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct vector {

  int dim;
  double *coeff;

} vector;

typedef struct state_vector {

  double x;
  double y;
  double z;
  double px;
  double py;
  double pz;

} state_vector;

typedef struct track {

  int num_steps;
  state_vector *states;

} track;

struct vector make_vector(int dim, double *coeff);

struct vector make_vector_2d(double c1, double c2);

struct vector make_vector_3d(double c1, double c2, double c3);

void scale_vector(double scale, struct vector *v);

void print_vector(const struct vector *v);

struct vector add_vectors(const struct vector *v1, const struct vector *v2);

struct vector add_scaled_vectors(double scale_1, const struct vector *v1, double scale_2, const struct vector *v2);

double scalar_product(const struct vector *v1, const struct vector *v2);

struct vector cross_product(const struct vector *v1, const struct vector *v2);

double compute_vector_norm(const struct vector *v);

struct state_vector add_scaled_state_vectors(double scale_1, const struct state_vector *v1, double scale_2, const struct state_vector *v2);

struct track make_empty_track(void);

void print_state(const struct state_vector *vec);

#endif