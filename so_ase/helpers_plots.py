import ipynbname
import matplotlib.pyplot as plt

def add_notebook_path_to_fig(fig, y_position=-0.05, fontsize=8):
    """
    Adds the absolute path of the Jupyter Notebook to an existing figure.

    Parameters:
    fig (matplotlib.figure.Figure): The existing figure to modify.
    y_position (float): The vertical position of the text (default: 0.01, near the bottom).
    fontsize (int): Font size of the path text (default: 8).
    """
    # Get absolute path of the Jupyter Notebook
    try:
        notebook_path = ipynbname.path()  # Get full notebook path
    except:
        notebook_path = "Notebook path not found"

    # Add notebook path as text at the lower end of the figure
    fig.text(0.5, y_position, f"Generated from: {notebook_path}", wrap=True, 
             horizontalalignment='center', fontsize=fontsize)
    

def remove_axes_frame(ax, left=True, right=True, top=True, bottom=True,
                      ticks=True, labels=True):
    """
    Remove or hide selected sides of a matplotlib Axes frame.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to modify.
    left, right, top, bottom : bool, optional
        If True, remove the corresponding spine.
    ticks : bool, optional
        If True, remove ticks on removed spines.
    labels : bool, optional
        If True, remove tick labels on removed spines.
    """

    sides = {
        "left": left,
        "right": right,
        "top": top,
        "bottom": bottom,
    }

    for side, remove in sides.items():
        if remove:
            ax.spines[side].set_visible(False)

    if ticks:
        if left:
            ax.yaxis.set_ticks_position("none")
        if bottom:
            ax.xaxis.set_ticks_position("none")

    if labels:
        if left:
            ax.set_yticklabels([])
        if bottom:
            ax.set_xticklabels([])

    return