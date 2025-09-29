import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from scipy import stats

mpl.use('Agg')
hist_filled = {'alpha': 0.3, 'histtype': 'stepfilled'}
hist_step = {'histtype': 'step'}
my_dpi = 128
fontsize = 24
fontsize_legend = 16
RED = "#D55E00"
BLUE = "#0072B2"
GREEN = "#009E73"
ORANGE = "#E69F00"
PURPLE = "#CC79A7"
YELLOW = "#F0E442"
CYAN = "#56B4E9"
WHITE = "#FFFFFF"
BLACK = "#000000"


def tex_f(x):
    if x == 0:
        s = "0.0"
    elif 0.001 < abs(x) < 10:
        s = f"{x:6.3f}"
    elif 10 <= abs(x) < 10000:
        s = f"{x:6.1f}"
    else:
        s = f"{x:6.2g}"
        if "e" in s:
            mantissa, exp = s.split('e')
            s = mantissa + '\\times 10^{' + str(int(exp)) + '}'
    return s


def plot_cat_box(ax, x_label, y_label, df, annotate=True, scale=""):
    filter_nan = np.isfinite(df[y_label])
    df_input = df[filter_nan].copy()
    sorted_x = sorted(df_input[x_label].unique())
    # One boxplot per category
    if scale == "log":
        df_input = df_input[df_input[y_label] > 0.0]
    sns.boxplot(data=df_input, x=x_label, y=y_label, order=sorted_x, ax=ax, notch=True,
                medianprops={'color': 'black'}, showfliers=False, width=0.9, log_scale=(scale == "log"))

    # Add the p-values for each pair of categories
    if annotate:
        x = {cat: group[y_label].values for cat, group in df_input.groupby(x_label, observed=True)}
        pairs = []
        for i in range(1, len(sorted_x)):
            pairs.append((sorted_x[i - 1], sorted_x[i]))
        if len(pairs) > 0:
            pvalues = []
            for pair in pairs:
                pvalues.append(stats.ttest_ind(x[pair[0]], x[pair[1]], alternative="two-sided").pvalue)
            formatted_pvalues = [f'$p={tex_f(pvalue)}$' for pvalue in pvalues]
            annotator = Annotator(ax, pairs, data=df_input, x=x_label, y=y_label, order=sorted_x)
            annotator.configure(loc='inside', verbose=0)
            annotator.set_custom_annotations(formatted_pvalues)
            annotator.annotate()
    ax.set_xlim(-0.5, len(sorted_x) - 0.5)
    ax.set_xticks(range(len(sorted_x)))
    # Count the number of samples in each category
    counts = df_input[x_label].value_counts().reindex(sorted_x, fill_value=0)
    ax.set_xticklabels([cat + f"\nn={counts[cat]}" for cat in sorted_x])
    if len(sorted_x) < 6:
        plt.xticks(rotation=0, ha="center")
    else:
        plt.xticks(rotation=45, ha="right")
    ax.margins(x=0.05)
    ax.set_xlabel("")


def plot_bar(ax, x_label, y_label, df, q=30):
    filter_nan = (np.isfinite(df[x_label]) & np.isfinite(df[y_label]))
    df_input = df[filter_nan].copy()
    df_input["xqcut"] = pd.qcut(df_input[x_label], q=min(q, len(df) // 2), duplicates="drop")

    df_gb = df_input.groupby("xqcut", observed=True).agg({x_label: "mean", y_label: "mean"}).reset_index()
    if len(df_gb) < 2:
        return
    ax.bar(np.arange(len(df_gb)),
           df_gb[y_label] - df_gb[y_label].min() + 0.05 * (df_gb[y_label].max() - df_gb[y_label].min()))


def plot_scatter_bins(ax, x_label, y_label, df, q=30):
    filter_nan = (np.isfinite(df[x_label]) & np.isfinite(df[y_label]))
    df_input = df[filter_nan].copy()
    df_input["xqcut"] = pd.qcut(df_input[x_label], q=min(q, len(df) // 2), duplicates="drop")
    gb = df_input.groupby("xqcut", observed=True)
    df_gb = gb.agg({x_label: "mean", y_label: "mean"}).reset_index()
    # Calculate standard deviation for error bars
    df_gb[f"{y_label}_std"] = gb[y_label].std().values
    if len(df_gb) < 10:
        return
    corr = df_gb[x_label].corr(df_gb[y_label])
    ax.plot(df_gb[x_label], df_gb[y_label], "o", alpha=1.0, label=f"r={corr:.2g}")
    ax.errorbar(df_gb[x_label], df_gb[y_label], yerr=df_gb[f"{y_label}_std"],
                fmt='none', alpha=0.5, ecolor="black", capsize=3)
    # Plot the linear regression line
    m, b = np.polyfit(df_gb[x_label], df_gb[y_label], 1)
    x = np.linspace(df_gb[x_label].min(), df_gb[x_label].max(), 100)
    y = m * x + b
    ax.plot(x, y, color="red", alpha=0.5, label=f"y={m:.2g}x+{b:.2g} (r={corr:.2g})")
    ax.legend()


def plot_bins_box(ax, x_label, y_label, df, q=30, scale=""):
    filter_nan = (np.isfinite(df[x_label]) & np.isfinite(df[y_label]))
    df_input = df[filter_nan].copy()
    qeff = min(q, len(df) // 2)
    df_input["xqcut"] = pd.qcut(df_input[x_label], q=qeff, duplicates="drop", labels=[f"Q{i}" for i in range(qeff)])
    plot_cat_box(ax, x_label="xqcut", y_label=y_label, df=df_input, annotate=False, scale=scale)


def plot_2d_histogram(ax, x, y, bins=30, cmap="Blues"):
    """
    Plot a 2D histogram on the given axes as contour plot.
    """
    # Clip to avoid extreme values, use the 1% and 99% quantiles
    x = np.clip(x, np.quantile(x, 0.01), np.quantile(x, 0.99))
    y = np.clip(y, np.quantile(y, 0.01), np.quantile(y, 0.99))
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins)
    X, Y = np.meshgrid(xedges[:-1], yedges[:-1], indexing='ij')
    # Gaussian smoothing
    hist = stats.gaussian_kde(np.vstack([x, y]), bw_method='scott')(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    # Use pcolormesh for better performance with large datasets
    ax.pcolormesh(X, Y, hist.T, cmap=cmap, shading='auto')
    # Add contour lines
    ax.contour(X, Y, hist.T, levels=10, colors='black', linewidths=0.5, linestyles='dashed')


def annotate_category(df, clade_1, clade_2, v="ratio"):
    t_1_max = df[f"{v}_{clade_1}"].quantile(0.75)
    t_1_min = df[f"{v}_{clade_1}"].quantile(0.25)
    t_2_max = df[f"{v}_{clade_2}"].quantile(0.75)
    t_2_min = df[f"{v}_{clade_2}"].quantile(0.25)

    def annotate_row(row):
        # Split the dataframe into four categories based on the ratio:
        if row[f"{v}_{clade_1}"] > t_1_max and row[f"{v}_{clade_2}"] > t_2_max:
            return f"Div in both"
        elif row[f"{v}_{clade_1}"] <= t_1_min and row[f"{v}_{clade_2}"] <= t_2_min:
            return "Stab"
        elif row[f"{v}_{clade_1}"] > t_1_max and row[f"{v}_{clade_2}"] <= t_2_min:
            return f"Div in {clade_1}"
        elif row[f"{v}_{clade_1}"] <= t_1_min and row[f"{v}_{clade_2}"] > t_2_max:
            return f"Div in {clade_2}"
        else:
            return "Other"

    return df.apply(lambda row: annotate_row(row), axis=1)
