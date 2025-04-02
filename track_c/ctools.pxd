cdef extern from "magnet2D.h":

	ctypedef struct magnet_2d:
		pass

	ctypedef struct configuration:
		pass

	magnet_2d make_magnet_2d(int num_core,
							const double *core_x,
							const double *core_y,
							const double *core_z,
                            int num_wing,
							const double *wing_x,
							const double *wing_y,
							const double *wing_z,
                            double B_goal,
							int type)

	configuration make_configuration(int num_magnets, const int *num_core, const int *num_wing,
                                        const double *x_core, const double *y_core, const double *z_core,
                                        const double *x_wing, const double *y_wing, const double *z_wing,
                                        const double *B_goal, const int *types)

cdef extern from "rktools.h":
	
	ctypedef struct particle:
		pass

	ctypedef struct track_limits:
		pass

	particle *make_particles(int num_particles, const double *m, const double *q, const double *y0)

	track_limits make_track_limits(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)

	void track_particles(double *ret_data, int num_particles, double h, double T_max, int subsampling, const track_limits *limits, const configuration *config, const particle *prtcls, int *num_samples)