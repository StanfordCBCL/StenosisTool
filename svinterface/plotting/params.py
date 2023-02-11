import matplotlib.pyplot as plt
import matplotlib as mpl


def set_params(font = "Times New Roman", size = 18, linewidth = 1):
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
    