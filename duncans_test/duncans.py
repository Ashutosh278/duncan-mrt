import pandas as pd
import numpy as np
from itertools import combinations
from scipy.stats import studentized_range


class DuncanResults:
    """
    A custom object to hold and format results from Duncan's Multiple Range Test.
    Display can be in mean-sorted or sequential order based on sorting parameter.
    """

    def __init__(self, header, critical_ranges, grouping, means, pvals, diffs, sorting='mean'):
        self.header = header
        self.critical_ranges = critical_ranges
        self.grouping = grouping
        self.means = means
        self.pvals = pvals
        self.diffs = diffs
        self.sorting = sorting  # Store sorting preference

    def _format_header(self):
        """Format the header section."""
        lines = []
        lines.append("=" * 70)
        lines.append(f"{'Duncan\'s Multiple Range Test Results':^70}")
        lines.append("=" * 70)
        lines.append("Test Parameters:")
        lines.append(f" α = {self.header['alpha']:.3f} | dfₑ = {self.header['df_error']:.0f} | "
                     f"MSₑ = {self.header['ms_error']:.3f} | CV = {self.header['cv']:.2f}%")
        lines.append("")
        return lines

    def _format_critical_ranges(self):
        """Format the critical ranges table."""
        lines = []
        lines.append("Critical Ranges:")
        lines.append("-" * 26)
        lines.append(f"{'Step':<6} {'Tprob':<8} {'Duncan':<8}")
        lines.append("-" * 26)
        for _, row in self.critical_ranges.round(3).iterrows():
            lines.append(f"{int(row['step']):<6} {row['Tprob']:<8.3f} {row['DUNCAN']:<8.3f}")
        lines.append("")
        return lines

    def _format_treatment_groups(self):
        lines = []
        # Apply sorting preference
        if self.sorting == 'mean':
            grouping_to_display = self.grouping.sort_values(by="Mean", ascending=False)
        else:
            grouping_to_display = self.grouping  # Sequential order

        # Dynamic column width calculation
        treatment_width = max(len("Treatment"),
                              grouping_to_display.index.astype(str).str.len().max() if not grouping_to_display.empty else 0)
        mean_width = max(len("Mean"),
                         len(f"{grouping_to_display['Mean'].max():.3f}") if not grouping_to_display.empty else 0)
        group_width = max(len("Group"),
                          grouping_to_display['Group'].astype(str).str.len().max() if not grouping_to_display.empty else 0)
        
        # Table formatting
        header_line = f"{'Treatment':<{treatment_width}}  {'Mean':<{mean_width}}  {'Group':<{group_width}}"
        separator = "-" * len(header_line)
        lines.append("Treatment Groups:")
        lines.append(separator)
        lines.append(header_line)
        lines.append(separator)
        
        # Data rows
        for index, row in grouping_to_display.iterrows():
            lines.append(f"{index:<{treatment_width}}  {row['Mean']:<{mean_width}.3f}  {row['Group']:<{group_width}}")
        lines.append("")
        lines.append("Treatments with the same letter are not significantly different")
        lines.append("=" * len(header_line))
        return lines

    def _repr_html_(self):
        """Generate HTML representation for Jupyter/Colab with superscripts."""
        # Start HTML structure
        html = """
        <div style="font-family: monospace; font-size: 12px;">
        <h3 style="text-align: center;">Duncan's Multiple Range Test Results</h3>
        <p><strong>Test Parameters:</strong><br>
        α = {alpha:.3f} | dfₑ = {df_error:.0f} | MSₑ = {ms_error:.3f} | CV = {cv:.2f}%
        </p>
        """.format(**self.header)

        # Critical Ranges Table
        html += """
        <p><strong>Critical Ranges:</strong></p>
        <table border="1" class="dataframe" style="border-collapse: collapse;">
          <thead>
            <tr style="text-align: right;">
              <th>Step</th>
              <th>Tprob</th>
              <th>Duncan</th>
            </tr>
          </thead>
          <tbody>
        """
        for _, row in self.critical_ranges.round(3).iterrows():
            html += f"<tr><td>{int(row['step'])}</td><td>{row['Tprob']:.3f}</td><td>{row['DUNCAN']:.3f}</td></tr>"
        html += """
          </tbody>
        </table>
        <br>
        """

        # Apply sorting preference for treatment groups
        if self.sorting == 'mean':
            display_grouping = self.grouping.sort_values(by="Mean", ascending=False)
        else:
            display_grouping = self.grouping  # Sequential order

        # Treatment Groups Table with Superscripts
        html += """
        <p><strong>Treatment Groups:</strong></p>
        <table border="1" class="dataframe" style="border-collapse: collapse;">
          <thead>
            <tr style="text-align: right;">
              <th>Treatment</th>
              <th>Mean</th>
              <th>Group</th>
            </tr>
          </thead>
          <tbody>
        """
        # Iterate in the chosen order
        for index, row in display_grouping.iterrows():
            mean_val = row['Mean']
            group_letters = row['Group']
            mean_sup_html = f"{mean_val:.3f}<sup>{group_letters}</sup>"
            html += f"<tr><td>{index}</td><td>{mean_sup_html}</td><td>{group_letters}</td></tr>"

        html += """
          </tbody>
        </table>
        <p><em>Treatments with the same letter are not significantly different.</em></p>
        </div>
        """
        return html

    def __str__(self):
        """
        Returns a formatted, human-readable string representation of the results.
        This method is automatically called by the print() function.
        """
        lines = []
        lines.extend(self._format_header())
        lines.extend(self._format_critical_ranges())
        lines.extend(self._format_treatment_groups())
        return "\n".join(lines)


def _order_pvalue(treatment, means, alpha, pvalue):
    """
    Orders treatments by mean and assigns grouping letters based on p-values.
    """
    n = len(means)
    # Expanded character set for groups
    letter = list("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") + [".", "+", "-", "*", "/", "#", "$", "%", "&", "^", "[", "]", ":", "@", ";", "_", "?", "!", "=", "#"]
    letter += [" "] * 2000 # Pad with spaces if needed
    df = pd.DataFrame({"treatment": treatment, "means": means})
    w = df.sort_values("means", ascending=False).reset_index()
    M = [""] * n
    k = 0
    j = 0
    new_group_flag = 0
    loop_counter = 0
    M[0] = letter[k]
    q = w["index"].to_numpy()
    while j < n - 1:
        loop_counter += 1
        if loop_counter > n:
            break
        for i in range(j, n):
            is_not_significant = pvalue[q[i], q[j]] > alpha
            if is_not_significant:
                if not M[i].endswith(letter[k]):
                    M[i] += letter[k]
            else:
                k += 1
                group_change_index = i
                new_group_flag = 0
                ja = j
                for jj in range(group_change_index, n):
                    M[jj] += ""
                M[group_change_index] += letter[k]
                for v in range(ja, group_change_index + 1):
                    if pvalue[q[v], q[group_change_index]] <= alpha:
                        j += 1
                        new_group_flag = 1
                    else:
                        break
                break
        if new_group_flag == 0:
            j += 1
    out = pd.DataFrame({
        "means": w["means"].to_numpy(),
        "groups": M
    }, index=w["treatment"])
    return out


def _weighted_mean(means, counts):
    """Calculates the weighted mean."""
    means = np.asarray(means, dtype=float)
    counts = np.asarray(counts, dtype=float)
    return np.sum(means * counts) / np.sum(counts)


def _compute_cv(ms_error, overall_mean):
    """Computes the Coefficient of Variation."""
    return (np.sqrt(ms_error) * 100.0) / overall_mean


def _tprob_duncan(ntr, df_error, alpha):
    """
    Calculates the critical values (T-probabilities) for Duncan's test.
    """
    Tprob = []
    for i in range(2, ntr + 1):
        Tprob.append(studentized_range.ppf((1 - alpha) ** (i - 1), i, df_error))
    return np.array(Tprob)


def duncan_test(data, groups, *, anova_table=None, n_rep=None, ms_error=None,
                df_error=None, alpha=0.05, sorting='mean'):
    """
    Performs Duncan's Multiple Range Test for comparing treatment means.

    Args:
        data (list, array): The response variable data (means if n_rep provided).
        groups (list, array): The categorical treatment groups.
        anova_table (pd.DataFrame, optional): ANOVA table.
        n_rep (int or list, optional): Number of replicates per mean.
        ms_error (float, optional): Mean Square Error.
        df_error (int, optional): Degrees of Freedom for Error.
        alpha (float, optional): Significance level. Defaults to 0.05.
        sorting (str): 'mean' for descending mean order, 'sequential' for input order.
                       Defaults to 'mean'.

    Returns:
        DuncanResults: A custom object containing the test results.
    """
    groups = np.asarray(groups)
    if anova_table is not None:
        if (ms_error is not None) or (df_error is not None):
            raise ValueError("Provide either anova_table OR (ms_error & df_error), not both.")
        anova = anova_table.copy()
        if "mean_sq" not in anova.columns and {"sum_sq", "df"}.issubset(anova.columns):
            anova["mean_sq"] = anova["sum_sq"] / anova["df"]
        for resid_name in ("Residual", "Residuals", "Error"):
            if resid_name in anova.index:
                ms_error = float(anova.loc[resid_name, "mean_sq"])
                df_error = float(anova.loc[resid_name, "df"])
                break
        else:
            ms_error = float(anova.iloc[-1]["mean_sq"])
            df_error = float(anova.iloc[-1]["df"])
        df = pd.DataFrame({"value": np.asarray(data, dtype=float), "treatment": groups})
        g = df.groupby("treatment", observed=True)
        means = g["value"].mean().values
        counts = g["value"].count().values
        treatments = g["value"].mean().index.values
        overall_mean = df["value"].mean()
    elif (n_rep is not None) and (ms_error is not None) and (df_error is not None):
        # Case where means and n_rep are provided directly
        means = np.asarray(data, dtype=float)
        treatments = groups
        if np.isscalar(n_rep):
            counts = np.repeat(int(n_rep), len(means)).astype(int)
        else:
            counts = np.asarray(n_rep, dtype=int)
        if len(counts) != len(means):
            raise ValueError("Length of n_rep must match number of means")
        overall_mean = _weighted_mean(means, counts)
    else:
        raise ValueError("Must provide either anova_table OR (n_rep, ms_error, df_error).")

    ntr = len(means)
    if len(treatments) != ntr:
        raise ValueError("Length mismatch: treatments and means.")

    unique_n = np.unique(counts)
    if len(unique_n) == 1:
        se_diff = np.sqrt(ms_error / unique_n[0])
    else:
        n_harm = ntr / np.sum(1.0 / counts)
        se_diff = np.sqrt(ms_error / n_harm)

    Tprob = _tprob_duncan(ntr, df_error, alpha)
    DUNCAN = Tprob * se_diff
    crit_df = pd.DataFrame({"step": np.arange(2, ntr + 1), "Tprob": Tprob, "DUNCAN": DUNCAN})

    # Get the order based on means for calculating p-values correctly
    order_idx = np.argsort(means)[::-1] # Descending order of means
    ord_index = np.argsort(order_idx)   # Original index positions

    pvals = np.ones((ntr, ntr), dtype=float)
    diffs = pd.DataFrame(index=treatments, columns=treatments, dtype=float)

    for i, j in combinations(range(ntr), 2):
        diff = abs(means[i] - means[j])
        # odif is the range in ranking positions + 1
        odif = abs(ord_index[i] - ord_index[j]) + 1
        stat = diff / se_diff
        # Calculate raw p-value using studentized range CDF
        p_raw = studentized_range.cdf(stat, odif, df_error)
        # Adjust p-value for Duncan's test
        p_adj = 1.0 - p_raw ** (1.0 / (odif - 1.0)) if odif > 1 else 1.0
        pvals[i, j] = pvals[j, i] = p_adj

    # Calculate mean differences matrix
    for a in range(ntr):
        for b in range(ntr):
            diffs.iloc[a, b] = means[a] - means[b]
    diffs.index = treatments
    diffs.columns = treatments

    # Assign group letters using the internal sorting logic (based on means)
    cld = _order_pvalue(treatments, means, alpha, pvals)

    # Create grouping DataFrame based on sorting preference
    if sorting == 'mean':
        grouping = cld[['means', 'groups']].rename(columns={'means': 'Mean', 'groups': 'Group'})
    else:  # 'sequential'
        gdf = pd.DataFrame({"Treatment": cld.index, "Mean": cld["means"].values, "Group": cld["groups"].values})
        grouping = gdf.set_index("Treatment").reindex(treatments)[["Mean", "Group"]]

    means_df = pd.DataFrame({"treatment": treatments, "mean": means, "n": counts})
    cv = _compute_cv(ms_error, overall_mean)

    header = {
        "alpha": alpha,
        "df_error": float(df_error),
        "ms_error": float(ms_error),
        "cv": float(cv)
    }

    return DuncanResults(header, crit_df, grouping, means_df, pvals, diffs, sorting=sorting)
