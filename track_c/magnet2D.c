#include "magnet2D.h"

struct magnet_2d make_magnet_2d(int num_core, const double *core_x, const double *core_y, const double *core_z,
                              int num_wing, const double *wing_x, const double *wing_y, const double *wing_z,
                              double B_goal, int type){

  // make the return object
  struct magnet_2d mag;

  // the core and the wing polygon
  mag.core = make_empty_polygon();
  mag.wing = make_empty_polygon();

  // fill them
  for (int i = 0; i < num_core; ++i){
    add_corner_point(&mag.core, core_x[i], core_y[i], core_z[i]);
    add_corner_point(&mag.wing, wing_x[i], wing_y[i], wing_z[i]);
  }

  // compute the areas
  double S_core = compute_polygon_surface(&mag.core);
  double S_wing = compute_polygon_surface(&mag.wing);

  if (type == 1){
    mag.B_core = B_goal;
    mag.B_wing = -1.0*S_core*B_goal/S_wing;
  }
  else if (type == 3){
    mag.B_wing = B_goal;
    mag.B_core = -1.0*S_wing*B_goal/S_core;
  }
  else{
    printf("magnet type %d unknown!\n", type);
  }

  return mag;

}

struct configuration make_empty_configuration(void){

  struct configuration config;

  config.num_magnets = 0;

  config.magnets = (struct magnet_2d*) calloc(0, sizeof(struct magnet_2d));

  return config;

}

void add_magnet_to_configuration(struct magnet_2d *mag, struct configuration *config){

  // the number of magnets
  int num_magnets = config->num_magnets + 1;

  // allocate the new space for the magnets
  struct magnet_2d *magnets_new = (struct magnet_2d*) calloc(num_magnets, sizeof(struct magnet_2d));

  for(int i = 0; i < num_magnets-1; ++i){

    magnets_new[i].core = config->magnets[i].core;
    magnets_new[i].wing = config->magnets[i].wing;
    magnets_new[i].B_core = config->magnets[i].B_core;
    magnets_new[i].B_wing = config->magnets[i].B_wing;
  }


  magnets_new[num_magnets-1].core = mag->core;
  magnets_new[num_magnets-1].wing = mag->wing;
  magnets_new[num_magnets-1].B_core = mag->B_core;
  magnets_new[num_magnets-1].B_wing = mag->B_wing;


  config->num_magnets = num_magnets;
  free(config->magnets);
  config->magnets = magnets_new;

  return;
}

int is_inside_magnet(const struct state_vector *state, const struct magnet_2d *mag){

  int check;

  check = is_inside_polygon(state, &mag->core);

  if (check == 1){
    return 1;
  }

  check = is_inside_polygon(state, &mag->wing);

  if (check == 1){
    return 2;
  }

  return 0;

}

// struct vector compute_B_from_configuration(const struct vector *point, const struct configuration *config){

//   int check;

//   for(int i = 0; i < config->num_magnets; ++i){
//     check = is_inside_magnet(point, &config->magnets[i]);

//     if(check == 1){
//       return compute_B(point, &config->magnets[i]);
//     }
//   }
//   return make_vector_3d(0.0, 0.0, 0.0);
// }

struct configuration make_configuration(int num_magnets, const int *num_core, const int *num_wing,
                                        const double *x_core, const double *y_core, const double *z_core,
                                        const double *x_wing, const double *y_wing, const double *z_wing,
                                        const double *B_goal, const int *types){

  // auxiliary variables
  int cnt_core = 0;  
  int cnt_wing = 0;                    
                                      
  // the configuration
  struct configuration config = make_empty_configuration();

  // allocate the new space for the magnets
  struct magnet_2d *magnets = (struct magnet_2d*) calloc(num_magnets, sizeof(struct magnet_2d));

  for(int i = 0; i < num_magnets; ++i){

    magnets[i] = make_magnet_2d(num_core[i], &x_core[cnt_core], &y_core[cnt_core], &z_core[cnt_core],
                                num_wing[i], &x_wing[cnt_wing], &y_wing[cnt_wing], &z_wing[cnt_wing],
                               B_goal[i], types[i]);
    
    cnt_core += num_core[i];
    cnt_wing += num_wing[i];

  }


  config.num_magnets = num_magnets;
  free(config.magnets); 
  config.magnets = magnets;

  return config;                      
}
