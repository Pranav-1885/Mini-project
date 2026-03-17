import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv
import os

# Read DP matrix CSV
with open("matrix.csv", newline="") as f:
    reader = list(csv.reader(f))

col_headers = reader[0]
row_headers = [row[0] for row in reader[1:]]
data        = [row[1:] for row in reader[1:]]

# ── Read traceback path
path_cells = set()
if os.path.exists("traceback_path.csv"):
    with open("traceback_path.csv", 
    newline="") as f:
        reader_p = csv.DictReader(f)
        for row in reader_p:
            r = int(row["row"])    # convert to data[] index
            c = int(row["col"])    # convert to data[] index
            path_cells.add((r, c))

# ── Read alignment strings 
align1 = align2 = ""
score  = ""
if os.path.exists("alignment.txt"):
    with open("alignment.txt") as f:
        lines = f.read().splitlines()
        if len(lines) >= 3:
            align1, align2, score = lines[0], lines[1], lines[2]

# ── Build figure 
fig, ax = plt.subplots(figsize=(max(8, len(col_headers) * 0.7 + 1),
                                max(5, len(row_headers) * 0.55 + 2)))
ax.axis("off")

table = ax.table(
    cellText=data,
    colLabels=col_headers,
    rowLabels=row_headers,
    loc="center",
    cellLoc="center"
)

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.4)

# ── Colour cells
TRACEBACK_COLOR = "#5DBB63"   # green
HEADER_COLOR    = "#D6E4F0"   # soft blue for headers
DEFAULT_COLOR   = "white"

for (row_idx, col_idx), cell in table.get_celld().items():
    if row_idx == 0 or col_idx == -1:
        cell.set_facecolor(HEADER_COLOR)
        cell.set_text_props(fontweight="bold")
    else:
        data_r = row_idx - 1
        data_c = col_idx 
        if (data_r, data_c) in path_cells:
            cell.set_facecolor(TRACEBACK_COLOR)
            cell.set_text_props(fontweight="bold", color="white")
        else:
            cell.set_facecolor(DEFAULT_COLOR)

    cell.set_edgecolor("#AAAAAA")

# ── Title & alignment annotation
plt.title("Smith-Waterman DP Matrix", fontsize=14, fontweight="bold", pad=12)

if align1 and align2:
    # Build a match line (| for match/partial, space for gap)
    match_line = ""
    for a, b in zip(align1, align2):
        if a == b:
            match_line += "|"
        elif a == "-" or b == "-":
            match_line += " "
        else:
            match_line += "·"   # mismatch / partial match

    alignment_text = (
        f"Optimal Local Alignment (score = {score})\n\n"
        f"  Seq 1 : {align1}\n"
        f"          {match_line}\n"
        f"  Seq 2 : {align2}"
    )

    fig.text(
        0.5, 0.02,
        alignment_text,
        ha="center", va="bottom",
        fontsize=10,
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#F0F8FF",
                  edgecolor="#AAAAAA", alpha=0.9)
    )

# ── Legend
legend_patch = mpatches.Patch(color=TRACEBACK_COLOR, label="Traceback path")
ax.legend(handles=[legend_patch], loc="upper right",
          bbox_to_anchor=(1.0, 1.08), fontsize=9, framealpha=0.8)

plt.tight_layout(rect=[0, 0.15, 1, 1])
plt.show()