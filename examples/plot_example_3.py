"""Example: Generate an annimation. This takes some time (few minutes).
   ==========
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime

from scipy.optimize import minimize

import trackship as ts

# %%
# Script parameters
# =================

# the log filename
LOG_FN = 'opt_log_2025-04-02_13-15-17.txt'

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

# %%
# Make the animation
# ==================
ts.animate_optimization('test.gif', LOG_FN, ts.type_1_generator, 6, gen_args)