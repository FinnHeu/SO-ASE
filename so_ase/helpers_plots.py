import ipynbname
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

###---> 1. General Appearance

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

###---> 2. Regression Lines and CIs
def plot_linear_trend_ci(ax, x, y, color='b', alpha=0.1, label='Linear trend'):
    """
    Plot linear regression with 95% confidence interval.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    x : array-like
        X data points.
    y : array-like
        Y data points.
    color : str, optional
        Color for the regression line and CI. Default is 'b'.
    alpha : float, optional
        Transparency of CI shading. Default is 0.1.
    label : str, optional
        Label for the regression line. Default is 'Linear trend'.
    
    Returns
    -------
    tuple
        slope, intercept, slope_ci : Regression parameters and confidence interval.
    """

    x_vals = np.asarray(x)
    y_vals = np.asarray(y)

    order = np.argsort(x_vals)
    x_vals = x_vals[order]
    y_vals = y_vals[order]

    n = len(x_vals)

    # Linear regression
    slope, intercept, r, p, stderr = stats.linregress(x_vals, y_vals)

    y_pred = intercept + slope * x_vals

    # Residuals
    residuals = y_vals - y_pred
    dof = n - 2

    # Residual std error
    s_err = np.sqrt(np.sum(residuals**2) / dof)

    # t-value for 95% CI
    t_val = stats.t.ppf(0.975, dof)

    x_mean = np.mean(x_vals)
    Sxx = np.sum((x_vals - x_mean)**2)

    ci = t_val * s_err * np.sqrt(
        1/n + (x_vals - x_mean)**2 / Sxx
    )

    ci_lower = y_pred - ci
    ci_upper = y_pred + ci

    # Plot regression
    label = f"{label}: {str(slope)[:6]} Mio. km2 per year"
    ax.plot(x_vals, y_pred, color=color, lw=.5, linestyle=':', label=label)

    ax.fill_between(
        x_vals,
        ci_lower,
        ci_upper,
        color=color,
        alpha=alpha,
        label='95% CI'
    )

    slope_ci = t_val * stderr

    return slope, intercept, slope_ci