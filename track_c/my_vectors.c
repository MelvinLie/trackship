#include "my_vectors.h"

struct vector make_vector(int dim, double *coeff){

  // allocate the return vector
  struct vector v_ret;
  
  // set the return dimension
  v_ret.dim = dim;

  // fill it
  v_ret.coeff = coeff;

  // return
  return v_ret;

}

void print_vector(const struct vector *v){

  printf("( ");
  for (int i = 0; i < v->dim-1; ++i){
    printf(" %f ,", v->coeff[i]);
  }
  printf(" %f )\n", v->coeff[v->dim-1]);

  return;

}

struct vector make_vector_2d(double c1, double c2){

  // allocate the return vector
  struct vector v_ret;

  // set the return dimension
  v_ret.dim = 2;

  // allocate the coefficients
  v_ret.coeff = (double*) calloc(2, sizeof(double));

  v_ret.coeff[0] = c1;
  v_ret.coeff[1] = c2;

  return v_ret;
  
}

struct vector make_vector_3d(double c1, double c2, double c3){

  // allocate the return vector
  struct vector v_ret;

  // set the return dimension
  v_ret.dim = 3;

  // allocate the coefficients
  v_ret.coeff = (double*) calloc(3, sizeof(double));

  v_ret.coeff[0] = c1;
  v_ret.coeff[1] = c2;
  v_ret.coeff[2] = c3;

  return v_ret;
}

void scale_vector(double scale, struct vector *v){

  // scale it
  for (int i = 0; i < v->dim; ++i){
    v->coeff[i] *= scale;
  }

}

struct vector add_vectors(const struct vector *v1, const struct vector *v2){

  // allocate the return vector
  struct vector v_ret;

  // check consistency
  if(v1->dim != v2->dim){
    printf("Vector dimensions are inconsistent!\n");
  }

  // set the return dimension
  v_ret.dim = v1->dim;

  // allocate the coefficients
  v_ret.coeff = (double*) calloc(v_ret.dim, sizeof(double));

  // fill it
  for (int i = 0; i < v1->dim; ++i){
    v_ret.coeff[i] = v1->coeff[i] + v2->coeff[i];
  }

  return v_ret;

}

struct vector add_scaled_vectors(double scale_1, const struct vector *v1, double scale_2, const struct vector *v2){

  // allocate the return vector
  struct vector v_ret;

  // check consistency
  if(v1->dim != v2->dim){
    printf("Vector dimensions are inconsistent!\n");
  }

  // set the return dimension
  v_ret.dim = v1->dim;

  // allocate the coefficients
  v_ret.coeff = (double*) calloc(v_ret.dim, sizeof(double));

  // fill it
  for (int i = 0; i < v1->dim; ++i){
    v_ret.coeff[i] = scale_1*v1->coeff[i] + scale_2*v2->coeff[i];
  }

  return v_ret;
}

double scalar_product(const struct vector *v1, const struct vector *v2){

  // the return value
  double ret_val = 0.0;

  // check consistency
  if(v1->dim != v2->dim){
    printf("Vector dimensions are inconsistent!\n");
  }

  // compute the scalar product
  for (int i = 0; i < v1->dim; ++i){
    ret_val += v1->coeff[i] * v2->coeff[i];
  }

  return ret_val;
}


struct vector cross_product(const struct vector *v1, const struct vector *v2){

  // allocate the return vector
  struct vector v_ret;

  // check consistency
  if(v1->dim != v2->dim){
    printf("Vector dimensions are inconsistent!\n");
  }
  
  if(v1->dim ==3){

    // set the return dimension
    v_ret.dim = 3;

    // allocate the coefficients
    v_ret.coeff = (double*) calloc(v_ret.dim, sizeof(double));

    v_ret.coeff[0] = v1->coeff[1]*v2->coeff[2] - v1->coeff[2]*v2->coeff[1];
    v_ret.coeff[1] = v1->coeff[2]*v2->coeff[0] - v1->coeff[0]*v2->coeff[2];
    v_ret.coeff[2] = v1->coeff[0]*v2->coeff[1] - v1->coeff[1]*v2->coeff[0];
    
  }
  else{

    printf("Cross product for dim != 3 not implemented!\n");
  }


  return v_ret;
}

double compute_vector_norm(const struct vector *v){

  // initialize the return value
  double norm = 0.0;

  for (int i = 0; i < v->dim; ++i){
    norm += v->coeff[i]*v->coeff[i];
  }

  return sqrt(norm);
}

struct state_vector add_scaled_state_vectors(double scale_1,
                                                     const struct state_vector *v1,
                                                     double scale_2,
                                                     const struct state_vector *v2){
  
  // allocate the return vector
  struct state_vector v_ret;
  
  v_ret.x  = scale_1*v1->x  + scale_2*v2->x;
  v_ret.y  = scale_1*v1->y  + scale_2*v2->y;
  v_ret.z  = scale_1*v1->z  + scale_2*v2->z;
  v_ret.px = scale_1*v1->px + scale_2*v2->px;
  v_ret.py = scale_1*v1->py + scale_2*v2->py;
  v_ret.pz = scale_1*v1->pz + scale_2*v2->pz;

  return v_ret;

}


struct track make_empty_track(void){

  // make the return value
  struct track tr;

  // allocate
  tr.num_steps = 0;
  tr.states =  (struct state_vector*) calloc(0, sizeof(struct state_vector));

  return tr;
}

void print_state(const struct state_vector *vec){

  printf("(x, y, z) = ( %e , %e, %e ) || (px, py, pz) = ( %f , %f, %f )\n", vec->x,
                                                                            vec->y,
                                                                            vec->z,
                                                                            vec->px/5.344286e-19,
                                                                            vec->py/5.344286e-19,
                                                                            vec->pz/5.344286e-19);

}