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