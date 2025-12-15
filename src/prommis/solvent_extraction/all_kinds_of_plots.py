from math import e
import pandas as pd
import matplotlib.pyplot as plt
from pyomo.environ import RangeSet

data_path = "feed_phi_parameters.xlsx"
df_alpha_1M_HCl = pd.read_excel(data_path, sheet_name="1 M HCl alpha values")
df_alpha_3M_HCl = pd.read_excel(data_path, sheet_name="3 M HCl alpha values")
df_alpha_5M_HCl = pd.read_excel(data_path, sheet_name="5 M HCl alpha values")

plt.plot(
    [df_alpha_1M_HCl.loc[i, "v"] for i in df_alpha_1M_HCl.index],
    [df_alpha_1M_HCl.loc[i, "alpha"] for i in df_alpha_1M_HCl.index],
    marker="o",
)
plt.plot(
    [df_alpha_3M_HCl.loc[i, "v"] for i in df_alpha_3M_HCl.index],
    [df_alpha_3M_HCl.loc[i, "alpha"] for i in df_alpha_3M_HCl.index],
    marker="o",
)
plt.plot(
    [df_alpha_5M_HCl.loc[i, "v"] for i in df_alpha_5M_HCl.index],
    [df_alpha_5M_HCl.loc[i, "alpha"] for i in df_alpha_5M_HCl.index],
    marker="o",
)
plt.xlabel("Feed volumetric flowrate (L/hr)")
plt.ylabel("Feed-membrane efficiency Pr")
plt.title("Pr feed-membrane efficiency profile")
plt.legend(["1 M HCl", "3 M HCl", "5 M HCl"])
plt.show()

df_alpha = pd.read_excel(data_path, sheet_name="All HCl data")
df_alpha = df_alpha.set_index(["Unnamed: 0"])

feed_values = [v for v in RangeSet(0.5, 10, 0.5)]
plt.plot(
    feed_values,
    [
        df_alpha.loc["Pr_c", "1 M HCl"]
        + df_alpha.loc["Pr_l", "1 M HCl"] * v
        + df_alpha.loc["Pr_q", "1 M HCl"] * v**2
        for v in feed_values
    ],
    marker="o",
)
plt.plot(
    feed_values,
    [
        df_alpha.loc["Pr_c", "3 M HCl"]
        + df_alpha.loc["Pr_l", "3 M HCl"] * v
        + df_alpha.loc["Pr_q", "3 M HCl"] * v**2
        for v in feed_values
    ],
    marker="o",
)
plt.plot(
    feed_values,
    [
        df_alpha.loc["Pr_c", "5 M HCl"]
        + df_alpha.loc["Pr_l", "5 M HCl"] * v
        + df_alpha.loc["Pr_q", "5 M HCl"] * v**2
        for v in feed_values
    ],
    marker="o",
)
plt.xlabel("Feed volumetric flowrate (L/hr)")
plt.ylabel("Feed-membrane efficiency Pr")
plt.title("Pr efficiency quadratic model profile")
plt.legend(["1 M HCl", "3 M HCl", "5 M HCl"])
plt.show()

df_alpha_alt = pd.read_excel(data_path, sheet_name="All HCl data alt")
df_alpha_alt = df_alpha_alt.set_index(["Unnamed: 0"])

HCl_conc = [1, 3, 5]

fig, ax = plt.subplots(3, 1, figsize=(5, 9))
fig.suptitle("Pr quadratic parameters variation")
ax[0].scatter(
    HCl_conc,
    [df_alpha_alt.loc["Pr_c", f"{e} M HCl"] for e in HCl_conc],
    marker="o",
)
ax[0].set(xlabel="HCl conc (M)", ylabel="Pr_c")
ax[0].set_xlim([0, 6])
ax[0].set_title("Pr_c")
ax[1].scatter(
    HCl_conc,
    [df_alpha_alt.loc["Pr_l", f"{e} M HCl"] for e in HCl_conc],
    marker="o",
)
ax[1].set(xlabel="HCl conc (M)", ylabel="Pr_l")
ax[1].set_xlim([0, 6])
ax[1].set_title("Pr_l")
ax[2].scatter(
    HCl_conc,
    [df_alpha_alt.loc["Pr_q", f"{e} M HCl"] for e in HCl_conc],
    marker="o",
)
ax[2].set(xlabel="HCl conc (M)", ylabel="Pr_q")
ax[2].set_title("Pr_q")
ax[2].set_xlim([0, 6])
plt.tight_layout()
