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
