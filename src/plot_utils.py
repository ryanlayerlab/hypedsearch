from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from matplotlib import rcParams

# Constants
SINGLE_FIG_SIZE = (6, 4)
BAR_WIDTH = 0.6
TICK_SIZE = 15
XLABEL_PAD = 10
LABEL_SIZE = 14
TITLE_SIZE = 16
LEGEND_SIZE = 11
LINE_WIDTH = 2
LIGHT_COLOR = "0.8"
LIGHT_COLOR_V = np.array([float(LIGHT_COLOR) for i in range(3)])
DARK_COLOR = "0.4"
DARK_COLOR_V = np.array([float(DARK_COLOR) for i in range(3)])
ALMOST_BLACK = "0.125"
ALMOST_BLACK_V = np.array([float(ALMOST_BLACK) for i in range(3)])
ACCENT_COLOR_1 = np.array([255.0, 145.0, 48.0]) / 255.0

# Configuration
# rcParams['text.usetex'] = True #Let TeX do the typsetting
rcParams["pdf.use14corefonts"] = True
rcParams["ps.useafm"] = True
# rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath'] #Force sans-serif math mode (for axes labels)
rcParams["font.family"] = "sans-serif"  # ... for regular text
# rcParams['font.sans-serif'] = ['Helvetica','Helvetica Neue', 'HelveticaNeue'] #, Avant Garde, Computer Modern Sans serif' # Choose a nice font here
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
rcParams["text.color"] = ALMOST_BLACK
rcParams["axes.unicode_minus"] = False

rcParams["xtick.major.pad"] = "8"
rcParams["axes.edgecolor"] = ALMOST_BLACK
rcParams["axes.labelcolor"] = ALMOST_BLACK
rcParams["lines.color"] = ALMOST_BLACK
rcParams["xtick.color"] = ALMOST_BLACK
rcParams["ytick.color"] = ALMOST_BLACK
rcParams["text.color"] = ALMOST_BLACK
rcParams["lines.solid_capstyle"] = "butt"


def single_fig(figsize=SINGLE_FIG_SIZE):
    return plt.subplots(1, 1, figsize=figsize)


def color_bp(bp, color):
    """Helper function for making prettier boxplots"""
    c = np.array(color)  # * 0.5
    c = tuple(c)

    for x in bp["boxes"]:
        plt.setp(x, color=c)
        x.set_facecolor(color)
    for x in bp["medians"]:
        plt.setp(x, color="w")
    for x in bp["whiskers"]:
        plt.setp(x, color=c)
    for x in bp["fliers"]:
        plt.setp(x, color=c)
    for x in bp["caps"]:
        plt.setp(x, color=c)


def adjust_spines(ax, spines):
    """From http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html"""
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(("outward", 10))  # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color("none")  # don't draw spine

    # turn off ticks where there is no spine
    if "left" in spines:
        ax.yaxis.set_ticks_position("left")
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if "bottom" in spines:
        ax.xaxis.set_ticks_position("bottom")
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


def hide_right_top_axis(ax):
    """Remove the top and right axis"""
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


# def set_font_size(font_size: int = LABEL_SIZE)


def finalize(axs, fontsize=LABEL_SIZE, labelpad=7, ignore_legend=False, add_grid=True):
    """Apply final adjustments"""
    try:
        ax = axs[0]
        for ax in axs:
            ax.tick_params(direction="out")
            hide_right_top_axis(ax)
            ax.yaxis.label.set_size(fontsize)
            ax.xaxis.label.set_size(fontsize)
            if ignore_legend == False:
                ax.legend(frameon=False)
            ax.tick_params(axis="both", which="major", labelsize=fontsize, pad=labelpad)
            if add_grid:
                ax.grid(True, linestyle="--", alpha=0.7)
    except:
        axs.tick_params(direction="out")
        hide_right_top_axis(axs)
        axs.yaxis.label.set_size(fontsize)
        axs.xaxis.label.set_size(fontsize)
        if ignore_legend == False:
            axs.legend(frameon=False)
        axs.tick_params(axis="both", which="major", labelsize=fontsize, pad=labelpad)


def lineswap_axis(fig, ax, zorder=-1000, lw=1, alpha=0.2, skip_zero=False):
    """Replace y-axis ticks with horizontal lines running through the background.
    Sometimes this looks really cool. Worth having in the bag 'o tricks.
    """
    fig.canvas.draw()  # Populate the tick vals/labels. Required for get_[xy]ticklabel calls.

    ylabels = [str(t.get_text()) for t in ax.get_yticklabels()]
    yticks = [t for t in ax.get_yticks()]
    [str(t.get_text()) for t in ax.get_xticklabels()]
    [t for t in ax.get_xticks()]

    x_draw = [
        tick for label, tick in zip(ylabels, yticks) if label != ""
    ]  # Which ones are real, though?
    y_draw = [tick for label, tick in zip(ylabels, yticks) if label != ""]

    xmin = x_draw[0]
    xmax = x_draw[-1]

    # Draw all the lines
    for val in y_draw:
        if val == 0 and skip_zero:
            continue  # Don't draw over the bottom axis
        ax.plot(
            [xmin, xmax],
            [val, val],
            color=ALMOST_BLACK,
            zorder=zorder,
            lw=lw,
            alpha=alpha,
        )

    ax.spines["left"].set_visible(False)  # Remove the spine
    ax.tick_params(axis="y", which="both", length=0)  # Erase ticks by setting length=0
    ax.set_xlim(xmin, xmax)  # Retain original xlims


def lighten_color(color, amount=0.5):
    """
        Lightens the given color by multiplying (1-luminosity) by the given amount.
        Input can be matplotlib color string, hex string, or RGB tuple.
    ​
        Examples:
        >> lighten_color('g', 0.3)
        >> lighten_color('#F034A3', 0.6)
        >> lighten_color((.3,.55,.1), 0.5)
    """
    import colorsys

    import matplotlib.colors as mc

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def fig_setup(nrows=1, ncols=1, w=6, h=4, constrained_layout: bool = True):
    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(ncols * w, nrows * h),
        constrained_layout=constrained_layout,
    )
    if nrows * ncols == 1:
        return fig, [axs]
    else:
        # return fig,axs.transpose().flatten()
        return fig, axs.flatten()


def plot_error_bars(ax, x, y, lower, upper, label=None, zorder=None, color=None):
    ax.errorbar(
        x,
        y,
        yerr=(y - lower, upper - y),
        fmt="o",
        label=label,
        zorder=zorder,
        color=color,
    )  # error is (lower, upper) = (mean-lower, upper-mean)
    return ax


def plot_centered_error_bars(ax, x, y, true_y, lower, upper, label=None, color=None):
    ax.errorbar(
        x, y - true_y, yerr=(y - lower, upper - y), fmt="o", label=label, color=color
    )  # error is (lower, upper) = (mean-lower, upper-mean)
    return ax


def set_title_axes_labels(ax, title=None, xlabel=None, ylabel=None):
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # return ax


def plot_pdf(ax, pdf, label=None, style="o"):
    ax.plot(pdf[:, 0], pdf[:, 1], style, label=label)


def plot_line(ax, m=1, b=0, label=None, ls="--", lc="black", lw=1):
    """
    Plot the y=m*x + b line
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    min_val = min(xlim[0], ylim[0])
    max_val = max(xlim[1], ylim[1])
    x = np.linspace(min_val, max_val, 100)
    y = m * x + b
    ax.plot(x, y, ls=ls, color=lc, lw=lw, label=label)


def best_fit_line(x, y):
    """
    Returns the slope and intercept of the line of best fit for the given x and y values.
    """
    # Calculate the slope and intercept using numpy's polyfit function
    slope, intercept = np.polyfit(x, y, 1)
    return slope, intercept


# Write me a function that (1) has a basic docstring that just describes what the function does
# in a few sentences, (2) computes the
def plot_best_fit_line(ax, x, y, label=None):
    """
    Plot the line of best fit for the given x and y values.
    """
    slope, intercept = best_fit_line(x, y)
    xlim = ax.get_xlim()
    ax.get_ylim()
    # min_val = min(xlim[0], ylim[0])
    # max_val = max(xlim[1], ylim[1])
    x = np.linspace(xlim[0], xlim[1], 1000)
    y = slope * x + intercept
    ax.plot(x, y, label=label)


def hist_plot(
    ax,
    data,
    nbins: int = 10,
    kde: bool = True,
    box_counts: bool = True,
    bw_adjust: float = 2,
    # label=None, color=None, alpha=0.5
):
    """
    Plot a histogram of the given data.
    """
    sns.histplot(
        data,
        bins=nbins,
        ax=ax,
        kde=kde,
        kde_kws=dict(bw_adjust=bw_adjust),
    )
    # if kde:
    #     sns.kdeplot(
    #         data,
    #         ax=ax,
    #         bw_adjust=2,
    #         label="KDE",
    #     )

    # Add count numbers to top of histogram boxes
    if box_counts:
        for patch in ax.patches:
            height = patch.get_height()
            if height > 0:  # only label non-empty bins
                _ = ax.text(
                    patch.get_x() + patch.get_width() / 2,
                    height,
                    f"{int(height)}",
                    ha="center",
                    va="bottom",
                )


def save_fig(
    path: Union[Path, str],
    dpi: int = 200,
    title: Optional[str] = None,
    fig: Optional[plt.Figure] = None,
    # bbox_inches: str = "tight",
    # transparent: bool = True,
    # fig_format: str = "pdf",
):
    """
    Save the figure to a file.
    """
    if title is not None:
        fig.suptitle(title, fontsize=16, fontweight="bold")
    plt.tight_layout()
    plt.savefig(
        path,
        dpi=dpi,
        bbox_inches="tight",
        # bbox_inches=bbox_inches,
        # transparent=transparent,
        # format=fig_format,
    )


@dataclass
class InteractiveScatterplot:
    df: pd.DataFrame
    x_colm: str
    y_colm: str
    title: str = ""
    color_colm: Optional[str] = None

    def plot(self, ms: int = 6, info_colms: Optional[List[str]] = None):
        """
        Create an interactive scatter plot using Plotly Express.
        """
        if info_colms is None:
            info_colms = list(self.df.columns)
        fig = px.scatter(
            self.df,
            x=self.x_colm,
            y=self.y_colm,
            color=self.color_colm,  # color by a column
            hover_data=info_colms,  # what to show on hover
        )
        _ = fig.update_traces(marker=dict(size=ms, opacity=0.7, line=dict(width=0)))
        fig.show()
        return fig


def interactive_scatter_plot(
    df: pd.DataFrame,
    x_colm: str,
    y_colm: str,
    title: str = "",
    color_colm: Optional[str] = None,
):
    info_colms = list(df.columns)
    fig = px.scatter(
        df,
        x=x_colm,
        y=y_colm,
        color=color_colm,  # color by a column
        hover_data=info_colms,  # what to show on hover
    )
    _ = fig.update_traces(
        marker=dict(size=6, opacity=0.7, line=dict(width=0))
    )  # improve visual clarity
    _ = fig.update_layout(
        title=title,
        hovermode="closest",
        xaxis=dict(title=x_colm),
        yaxis=dict(title=y_colm),
    )
    fig.show()
    return fig
