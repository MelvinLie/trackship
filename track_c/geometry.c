#include "geometry.h"


struct polygon make_empty_polygon(void){

  // make the return polygon
  struct polygon ret_poly;

  // set it up
  ret_poly.num_corners = 0;
  ret_poly.cp = (struct vector*) calloc(0, sizeof(struct vector));

  return ret_poly;

}

struct polygon make_polygon(int num_corners, const double *cp_x, const double *cp_y, const double *cp_z){

  // make the return polygon
  struct polygon ret_poly;

  // set it up
  ret_poly.num_corners = num_corners;
  ret_poly.cp = (struct vector*) calloc(num_corners, sizeof(struct vector));

  // fill the corner points
  for (int i = 0; i < num_corners; ++i){

    ret_poly.cp[i] = make_vector_3d(cp_x[i], cp_y[i], cp_z[i]);

  }
  
  return ret_poly;

}

void add_corner_point(struct polygon *poly, double x, double y, double z){

  // new number of corner points
  int num_corners = poly->num_corners + 1;

  // allocate the new space for the corner points
  struct vector *cp = (struct vector*) calloc(num_corners, sizeof(struct vector));

  // fill the corner points
  for (int i = 0; i < num_corners-1; ++i){
    cp[i] = make_vector_3d(poly->cp[i].coeff[0], poly->cp[i].coeff[1], poly->cp[i].coeff[2]);
  }

  // append the new point
  cp[num_corners-1] = make_vector_3d(x, y, z);

  // free the memory of the old data
  free(poly->cp);

  // set the new data
  poly->cp = cp;
  poly->num_corners = num_corners;

  return;

}

void print_polygon(const struct polygon *poly){

  printf("polygon with %d corners:\n", poly->num_corners);

  for(int i = 0; i < poly->num_corners; ++i){

    print_vector(&poly->cp[i]);
    
  }

}

double compute_polygon_surface(const struct polygon *poly){


  // the surface
  double S = 0.0;

  for (int i = 0; i < poly->num_corners-1; ++i){

    S += poly->cp[i].coeff[0]*poly->cp[i+1].coeff[2] - poly->cp[i+1].coeff[0]*poly->cp[i].coeff[2];

  }

  S += poly->cp[poly->num_corners-1].coeff[0]*poly->cp[0].coeff[2];
  S -= poly->cp[0].coeff[0]*poly->cp[poly->num_corners-1].coeff[2];

  S *= 0.5;

  return S;

}

int is_inside_polygon(const struct state_vector *state, const struct polygon *poly){

  // the solid angle
  double sa = 0.0;

  double d_1[3];
  double d_2[3];

  // the norms of the two vectors
  double norm_d_1, norm_d_2;

  // the scalar product of the two distance vectors
  double sp;

  // the argument of the arccos
  double arg;

  // loop over all corner points
  for (int i = 0; i < poly->num_corners; ++i){

    d_1[0] = poly->cp[i].coeff[0] - state->x;
    d_1[1] = poly->cp[i].coeff[1] - state->y;
    d_1[2] = poly->cp[i].coeff[2] - state->z;

    if (i < poly->num_corners - 1){
      d_2[0] = poly->cp[i + 1].coeff[0] - state->x;
      d_2[1] = poly->cp[i + 1].coeff[1] - state->y;
      d_2[2] = poly->cp[i + 1].coeff[2] - state->z;

      // d_2 = add_scaled_vectors(1.0, &poly->cp[i + 1], -1.0, point);
    }
    else{
      d_2[0] = poly->cp[0].coeff[0] - state->x;
      d_2[1] = poly->cp[0].coeff[1] - state->y;
      d_2[2] = poly->cp[0].coeff[2] - state->z;
      // d_2 = add_scaled_vectors(1.0, &poly->cp[0], -1.0, point);
    }

    norm_d_1 = sqrt(d_1[0]*d_1[0] + d_1[1]*d_1[1] + d_1[2]*d_1[2]);
    norm_d_2 = sqrt(d_2[0]*d_2[0] + d_2[1]*d_2[1] + d_2[2]*d_2[2]);
    // norm_d_1 = compute_vector_norm(&d_1);
    // norm_d_2 = compute_vector_norm(&d_2);

    // check already if we collapse on one of the points
    if (norm_d_1 < 1e-10){
      return 1;
    }
    if (norm_d_2 < 1e-10){
      return 1;
    }

    sp = d_1[0]*d_2[0] + d_1[1]*d_2[1] + d_1[2]*d_2[2];
    // sp = scalar_product(&d_1, &d_2);

    arg = sp/norm_d_1/norm_d_2;

    // numerical istabilities
    if (arg > 1.0){
      arg = 1.0;
    }
    else if (arg < -1.0){
      arg = -1.0;
    }

    sa += acos(arg);

  }

  if (fabs(sa - 2.0*M_PI) < 1e-6){

    return 1;
  }

  return 0;
}