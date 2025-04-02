import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

def animate_optimization(out_filename, log_filename, generator, num_var, args):
    '''Generate an animation showing all optimization steps, based on a log file.

    :param out_filename:
        The filename for the output.

    :param log_filename:
        The log filename to read.

    :param generator:
        The configuration generator.

    :param num_var:
        The number of variables.

    :param args:
        The additional generator arguments.
        
    :return:
        None
    '''

    # read the log file
    log_data = pd.read_csv(log_filename)

    # the number of optimization steps
    num_steps = log_data.values.shape[0]

    # make the configuration
    config = generator(log_data.values[0, :num_var], args)

    # the maximum values for the cost function
    max_fs = max(log_data.values[:, 6]**2)
    max_fc = max(log_data.values[:, 7])
    max_fa = max(log_data.values[:, 8])

    # the minimum values for the cost function
    min_fs = min(log_data.values[:, 6]**2)
    min_fc = min(log_data.values[:, 7])
    min_fa = min(log_data.values[:, 8])

    # the cost function values at the beginning
    f_s = log_data.values[0, 6]**2
    f_a = log_data.values[0, 7]
    f_c = log_data.values[0, 8]

    # open the figure
    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(221)
    ax1.set_xlim((-5, 5))
    ax1.set_ylim((-2, 38))
    ax1.set_title('footprint')

    ax2 = fig.add_subplot(222)
    ax2.bar(np.array([1, 2, 3])+0.2, [max_fs, max_fc, max_fa], width=0.2, color='red', label='max')
    ax2.bar(np.array([1, 2, 3])-0.2, [min_fs, min_fc, min_fa], width=0.2, color='green', label='min')
    bb = ax2.bar([1, 2, 3], [f_s, f_a, f_c], width=0.2, color='blue', label='current')
    ax2.set_xticks(np.array([1, 2, 3]), ['score', 'area', 'constraint'])
    ax2.set_yscale('log')
    ax2.legend()
    ax2.set_title('cost function')

    ax3 = fig.add_subplot(223)
    ax3.plot(np.linspace(1, num_steps, num_steps),
            log_data.values[:, -1], color='gray', linewidth=0.5)
    sc = ax3.plot(1, log_data.values[0, -1], 'o', color='red')[0]
    ax3.set_yscale('log')
    ax3.set_title('objective')

    # this is the animation body
    def animate(i):
        ax1.clear()
        sc.set_xdata(np.array([i+1]))
        sc.set_ydata(np.array([log_data.values[i, -1]]))
        config = generator(log_data.values[i, :num_var], args)
        config.plot(ax1)
        f_s = log_data.values[i, 6]**2
        f_a = log_data.values[i, 7]
        f_c = log_data.values[i, 8]
        bb.get_children()[0].set_height(f_s)
        bb.get_children()[1].set_height(f_a)
        bb.get_children()[2].set_height(f_c)
        ax1.set_xlim((-5, 5))
        ax1.set_ylim((-2, 40))
        ax1.set_title('footprint')
        ax2.set_title('cost function')
        ax3.set_title('objective')

        return

    ani = animation.FuncAnimation(fig, animate, repeat=False, frames=num_steps, interval=50)
    ani.save(out_filename, dpi=300, writer=PillowWriter(fps=25))
    plt.show()