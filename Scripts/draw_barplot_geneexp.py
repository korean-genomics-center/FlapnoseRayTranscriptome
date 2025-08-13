# %%
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

plt.rcParams["font.size"] = 16
# Select top 10 upregulated and top 10 downregulated based on log2fc
df_exp_cnt_filt["abslog2fc"] = df_exp_cnt_filt["log2fc"].abs()
df_plot = df_exp_cnt_filt[df_exp_cnt_filt["abslog2fc"] > 0.5]
top_up = df_plot.sort_values("log2fc", ascending=False).head(10)
top_down = df_plot.sort_values("log2fc", ascending=True).head(10)

# Set colors
top_up["Color"] = "#d73027"
top_down["Color"] = "#4575b4"

# Sort for plotting
top_up = top_up.sort_values("log2fc")
top_down = top_down.sort_values("log2fc")

# Create subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 10), facecolor="white")

# Add subplot labels A and B to the left side of the figure
fig.text(0.01, 0.95, "A", fontsize=22, fontweight='bold')  # Corresponds to top subplot
fig.text(0.01, 0.47, "B", fontsize=22, fontweight='bold')  # Corresponds to bottom subplot

# Plot upregulated (A)
axs[0].barh(top_up["gene_label"], top_up["log2fc"], color=top_up["Color"])
axs[0].set_title("Up-regulated", fontsize=18)
axs[0].set_xlabel("log2(Fold Change)")
axs[0].set_xlim(-3, 3)
axs[0].axvline(0, color="gray", linestyle="--")
axs[0].grid(axis='x', linestyle='--', alpha=0.5)

# Plot downregulated (B)
axs[1].barh(top_down["gene_label"], top_down["log2fc"], color=top_down["Color"])
axs[1].set_title("Down-regulated", fontsize=18)
axs[1].set_xlabel("log2(Fold Change)")
axs[1].set_xlim(-3, 3)
axs[1].axvline(0, color="gray", linestyle="--")
axs[1].grid(axis='x', linestyle='--', alpha=0.5)

# Formatting
for ax in axs:
    ax.set_xlim(-3, 3)
    ax.xaxis.set_major_locator(MultipleLocator(1))

plt.tight_layout()
plt.subplots_adjust(hspace=0.3)
plt.show()
plt.close()

# %%
