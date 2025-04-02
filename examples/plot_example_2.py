"""Example: Plot the optimization solution. Also generate the discretization.
   ==========
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime
from scipy import interpolate

from scipy.optimize import minimize

import trackship as ts

# %%
# Script parameters
# =================

# what is 1 GEV
GEV = 5.344286e-19

# the limit in the x direction for tracking
XTRACKLIM = 4.0

# the limit in the z direction for tracking
ZTRACKLIM = 35

# maximum number of steps
MAXSTEPS = 100000

# subsampling rate to reduce memory
SUBSR = 5

# the scoring plane width
WSP = 2.0

# arguments for the generator
gen_args = {
                "l_mhs(m)": 2.3,
                "w_mhs(m)": 0.3,
                "v_mhs(m)": 0.3,
                "b_mhs(T)": 1.9,
                "b_s1(T)": 1.8,
                "b_s2(T)": 1.7,
                "l_tot(m)": 29.7,
                "l_1(m)": 1.9*17.0/1.8,
                "dz(m)": 0.2
}


# the solution filename
solution_filename = 'opt_log_2025-04-02_13-15-17.npy'

# read the solution
# x = np.load(solution_filename)
x = np.array([0.115, 0.400, 0.100, 0.135, 0.300, 0.320])

# %%
# Make the particles
# ==================
NUM_P = 100
NUM_PHI = 9
PHI_MAX = 1*np.pi/180.0

NUM_MUONS = NUM_P*NUM_PHI

pz_0 = np.logspace(0.0, np.log10(350), NUM_P)*GEV
phi_0 = np.linspace(-PHI_MAX, PHI_MAX, NUM_PHI)

particles = []

all_phi = np.zeros((NUM_P*NUM_PHI, ))
all_p = np.zeros((NUM_P*NUM_PHI, ))

for i, pphi in enumerate(phi_0):
    for j, pp in enumerate(pz_0):

        particles.append(ts.Particle(105.66*1.79e-30,
                        1.0,
                        np.array([0., 0., 0., np.sin(pphi)*pp, 0., np.cos(pphi)*pp])))

        all_phi[i*NUM_P+j] = pphi
        all_p[i*NUM_P+j] = pp

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist2d(all_p/GEV, np.sin(all_phi)*all_p/GEV, norm=mpl.colors.LogNorm(), bins=100)
ax.set_xlabel('$|P|$ in GEV')
ax.set_ylabel('$P_t$ in GEV')
plt.show()

# %%
# Make the configuration
# =======================
config = ts.type_1_generator(x, gen_args)

# %%
# Make the discretization
# =======================

config_disc, _ = ts.type_1_discretizer(x, 4, 2, gen_args, plot=True)


# %%
# Track the muons
# ===============

# make the tracker for the original configuration
tracker = ts.RK4_ParticleTracker(config)

# make the tracker for the discretized configuration
tracker_disc = ts.RK4_ParticleTracker(config_disc)

# set the limits
tracker.set_x_limits(-XTRACKLIM, XTRACKLIM)
tracker.set_z_limits(-1.0, ZTRACKLIM)

tracker_disc.set_x_limits(-XTRACKLIM, XTRACKLIM)
tracker_disc.set_z_limits(-1.0, ZTRACKLIM)

# track
tracks = tracker.track(particles, MAXSTEPS, subsampling_rate=SUBSR, w_sp=WSP)
tracks_disc = tracker_disc.track(particles, MAXSTEPS, subsampling_rate=SUBSR, w_sp=WSP)

# %%
# Plot
# ===============
# number of colors
num_col = 30

# colormap for the muons according to their energy
colors = plt.cm.jet(np.linspace(0.0, 1.0, num_col))
f_r = interpolate.interp1d(np.logspace(-3, 0.0, num_col), colors[:, 0])
f_g = interpolate.interp1d(np.logspace(-3, 0.0, num_col), colors[:, 1])
f_b = interpolate.interp1d(np.logspace(-3, 0.0, num_col), colors[:, 2])

cmap = mpl.cm.ScalarMappable(cmap=mpl.cm.jet)
cmap.set_array([])

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
config.plot(ax1)
config_disc.plot(ax2)
for i in range(NUM_MUONS):

    color = np.array([f_r(all_p[i]/350/GEV), f_g(all_p[i]/350/GEV), f_b(all_p[i]/350/GEV), 0.1])

    tracks[i].plot(ax1, color=color)
    tracks_disc[i].plot(ax2, color=color)

cbar = plt.colorbar(cmap, ax=plt.gca(), ticks=[0.0, 0.333, 0.666, 1.0])
cbar.set_ticks(ticks=[0.0, 0.333, 0.666, 1.0])
cbar.set_ticklabels([f"{10**(-3)*350:.2f} GeV", f"{10**(-2)*350:.2f} GeV", f"{10**(-1)*350:.2f} GeV", f"{350:.2f} GeV"])

ax1.set_xlabel('$x$ in m')
ax1.set_ylabel('$y$ in m')
ax2.set_xlabel('$x$ in m')
ax2.set_ylabel('$y$ in m')
plt.show()

# %%
# Get snoopy dict
# ===============

# the vertical gap size
y_void = 0.5

# make an empty snoopy parameters dict
snoopy_params = {
    "yoke_type": [],
    "coil_material": [],
    "max_turns": [],
    "J_tar(A/mm2)": [],
    "NI(A)": [],
    "coil_diam(mm)": [],
    "insulation(mm)": [],
    "yoke_spacer(mm)": [],
    "material": [],
    "field_density": [],
    "delta_x(m)": [],
    "delta_y(m)": [],
    "delta_z(m)": [],
    "Z_pos(m)": [],
    "Z_len(m)": [],
    "Xmgap1(m)": [],
    "Xmgap2(m)": [],
    "Xcore1(m)": [],
    "Xcore2(m)": [],
    "Xvoid1(m)": [],
    "Xvoid2(m)": [],
    "Xyoke1(m)": [],
    "Xyoke2(m)": [],
    "Ycore1(m)": [],
    "Ycore2(m)": [],
    "Yvoid1(m)": [],
    "Yvoid2(m)": [],
    "Yyoke1(m)": [],
    "Yyoke2(m)": []
}


for i in range(len(config_disc.magnets)):
    if i % 2 == 0:
        config_disc.magnets[i].to_snoopy(snoopy_params, y_void)

snoopy_df = pd.DataFrame.from_dict(snoopy_params)

print(snoopy_df)

snoopy_df.to_csv('snoopy_parameters.csv', index=False)