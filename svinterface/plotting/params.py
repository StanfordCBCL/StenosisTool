import matplotlib.pyplot as plt
import matplotlib as mpl


def set_params(font = "serif", size = 12, linewidth = 1, use_latex = False, small_ticks = False):
    """Converts global parameters

    Args:
        font (str, optional): font. Defaults to "Helvetica Neue".
        size (int, optional): size of font. Defaults to 18.
        linewidth (int, optional): width of axes lines. Defaults to 2.
    """
    # Edit the font, font size, and axes width
    mpl.rcParams['font.family'] = font
    plt.rcParams['font.size'] = size
    plt.rcParams['axes.linewidth'] = linewidth
    if small_ticks:
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
    plt.rc('text', usetex=use_latex)