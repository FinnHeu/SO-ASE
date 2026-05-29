import numpy as np
import xarray as xr
import pandas as pd
from eofs.standard import Eof
from datetime import datetime
import matplotlib.pyplot as plt

def plot_ZW3_index(out_ds):
    """
    Plot Zonal Wave 3 (ZW3) index components with 3-month rolling means.
    
    Parameters
    ----------
    out_ds : xr.Dataset
        Output from ZW3_index(), must contain
        'pc1', 'pc2', 'zw3_magnitude', 'zw3_phase'.
    """

    # Extract PCs and indices
    pc1 = out_ds['pc1']
    pc2= out_ds['pc2']
    magnitude_index = out_ds['zw3_magnitude']
    phase_index = out_ds['zw3_phase']

    # Compute 3-month rolling means
    pc1_roll = pc1.rolling(time=3, center=True).mean()
    pc2_roll = pc2.rolling(time=3, center=True).mean()
    mag_roll = magnitude_index.rolling(time=3, center=True).mean()
    phase_roll = phase_index.rolling(time=3, center=True).mean()

    # Time vector
    time = pc1.time
    time_min, time_max = time_min, time_max = time[0], time[-1]

    # Function to add positive/negative shading
    def shade_pos_neg(ax, data, color_pos='red', color_neg='blue', alpha=0.8):
        ax.fill_between(time, 0, data, where=(data>0), facecolor=color_pos, alpha=alpha, interpolate=True)
        ax.fill_between(time, 0, data, where=(data<0), facecolor=color_neg, alpha=alpha, interpolate=True)

    # --- Plot setup ---
    fig, axes = plt.subplots(4, 1, figsize=(15, 10), sharex=True)

    # PC1
    axes[0].plot(time, pc1, color='grey', label='PC1', linewidth=0.5)
    axes[0].plot(time, pc1_roll, color='k', label='PC1 3-month rolling mean')
    shade_pos_neg(axes[0], pc1_roll)
    axes[0].set_ylabel('PC1 (normalized)')
    axes[0].set_ylim(-2.5, 2.5)
    axes[0].set_xlim(time_min, time_max)
    axes[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # PC2
    axes[1].plot(time, pc2, color='grey', label='PC2', linewidth=0.5)
    axes[1].plot(time, pc2_roll, color='k', label='PC2 3-month rolling mean')
    shade_pos_neg(axes[1], pc2_roll)
    axes[1].set_ylabel('PC2 (normalized)')
    axes[1].set_ylim(-2.5, 2.5)
    axes[1].set_xlim(time_min, time_max)
    axes[1].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # ZW3 magnitude
    axes[2].plot(time, magnitude_index, color='grey', label='ZW3 Magnitude', linewidth=0.5)
    axes[2].plot(time, mag_roll, color='k', label='ZW3 Magnitude 3-month rolling mean')
    axes[2].fill_between(time, mag_roll, 1.5, color='grey', alpha=0.95, where=(mag_roll > 1.5))
    axes[2].set_ylabel('ZW3 index (magnitude)')
    axes[2].set_xlim(time_min, time_max)
    axes[2].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # ZW3 phase
    axes[3].plot(time, phase_index, color='grey', label='ZW3 Phase', linewidth=0.5)
    axes[3].plot(time, phase_roll, color='k', label='ZW3 Phase 3-month rolling mean')
    axes[3].set_ylabel('ZW3 phase (deg)')
    axes[3].set_xlabel('Time')
    axes[3].set_ylim(-175, 175)
    axes[0].set_xlim(time_min, time_max)
    axes[3].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    plt.tight_layout()
    plt.show()
    
def plot_nino34_index(nino34_index, title="Niño 3.4 SST Anomaly Index", highlight=True, threshold=0.5):
    """
    Plot the Niño 3.4 index with ENSO event shading and labels.

    Parameters
    ----------
    nino34_index : xr.DataArray
        Output from nino34_index()
    title : str
        Plot title
    highlight : bool
        Whether to shade El Niño / La Niña periods
    threshold : float
        Threshold for moderate ENSO events (default 0.5°C)
    """
    
    # Convert to pandas for easier plotting
    ts = nino34_index.to_pandas()

    plt.figure(figsize=(14,4))
    plt.plot(ts.index, ts.values, label='Niño 3.4 SST Anomaly', color='black', linewidth=1.5)

    # Plot thresholds
    plt.axhline(threshold, linestyle='--', color='red', linewidth=1)
    plt.axhline(-threshold, linestyle='--', color='blue', linewidth=1)

    if highlight:
        # Moderate events shading
        plt.fill_between(ts.index, ts.values, threshold, where=(ts.values >= threshold), color='red', alpha=0.5)
        plt.fill_between(ts.index, ts.values, -threshold, where=(ts.values <= -threshold), color='blue', alpha=0.5)

    plt.title(title)
    plt.ylabel("SST Anomaly (°C)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.xlim(ts.index[0], ts.index[-1])
    plt.tight_layout()
    plt.show()