import numpy as np
import matplotlib.pyplot as plt

from .tracker import RK4_ParticleTracker

def obj_fcn_1(x, generator, particles, obj_args, config_args, plt_ax=None):
    '''The type 1 objective function for the muon shield optimization.
    The objective function is
        f = s**2 + a + c
    where
        s: Muon score,
        a: Total surface spanned by the magnets,
        c: Constraint (caver walls).
    
    :param x:
        The design variables:

    :param generator:
        The generator to use.

    :param particles:
        A list of particles to be used.

    :param obj_args:
        Additional arguments passed to the objective function.
        A dictionary with the following parameters:
        "x_track_min(m)" minimum x coordinate to continue tracking a particle.
        "x_track_max(m)" maximum x coordinate to continue tracking a particle.
        "z_track_min(m)" minimum z coordinate to continue tracking a particle.
        "z_track_max(m)" maximum z coordinate to continue tracking a particle.
        "max_steps" maximum steps to track.
        "subsampling_rate" the subsampling rate for memory reduction.
        "w_scoring_plane(m)" the width of the scoring plane.
        "x_lim(m)" the width of the cavern.
        "log_filename" the name of the log file to store the parameters.

    :param config_args:
        Additional arguments passed to the configuration generator.

    :param plt_ax:
        The plot axes. Default None. If no axes object is passed,
        no plot is generated.

    :return:
        The value of the objective function.
    '''

    # make the configuration
    config = generator(x, config_args)

    # make the tracker
    tracker = RK4_ParticleTracker(config)

    # set the limits
    tracker.set_x_limits(obj_args["x_track_min(m)"], obj_args["x_track_max(m)"])
    tracker.set_z_limits(obj_args["z_track_min(m)"], obj_args["z_track_max(m)"])

    # track
    tracks = tracker.track(particles, obj_args["max_steps"],
                           subsampling_rate=obj_args["subsampling_rate"],
                           w_sp=obj_args["w_scoring_plane(m)"])

    # count the hits
    s = 0
    for _, tr in enumerate(tracks):
        if tr.hit:
            s += abs(abs(tr.states[-1, 0]) - obj_args["w_scoring_plane(m)"])

    # the total surface
    a = config.compute_total_surface()

    # the constraint (cavern wall)
    c = 0.0

    # section 1
    if x[1]+2*x[0] > obj_args["x_lim(m)"]:
        c += (x[1]+2*x[0] - obj_args["x_lim(m)"])**2
    # section 2 (entrance)
    if 2*x[2]+x[3] > obj_args["x_lim(m)"]:
        c += (2*x[2]+x[3] - obj_args["x_lim(m)"])**2
    # section 2 (exit)
    if 2*x[4]+x[5] > obj_args["x_lim(m)"]:
        c += (2*x[4]+x[5] - obj_args["x_lim(m)"])**2

    with open(obj_args["log_filename"], "a") as myfile:
        myfile.write(f"{x[0]}")
        for i, xx in enumerate(x[1:]):
            myfile.write(f",{xx}")
        myfile.write(f",{s},{a},{c},{s*s+a+c}\n")

    print(f"s = {s:.2e}, a = {a:.2e}, c = {c:.2e}, f = {s*s + a + c:.2e}")

    if isinstance(plt_ax, plt.Axes):
        plt_ax.clear()
        config.plot(plt_ax)
        for i in range(len(particles)):
            tracks[i].plot(plt_ax)

        plt_ax.set_xlabel('$x$ in m')
        plt_ax.set_ylabel('$y$ in m')
        plt.pause(0.05)

    return s*s + a + c