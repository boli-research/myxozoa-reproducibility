from __future__ import annotations

from collections.abc import Mapping

from Bio.Phylo.BaseTree import Tree


def tree_layout(tree: Tree) -> tuple[dict, dict, dict[str, float]]:
    leaves = tree.get_terminals()
    leaf_y = {leaf: float(len(leaves) - 1 - index) for index, leaf in enumerate(leaves)}
    y = {}

    def assign_y(clade):
        if clade in leaf_y:
            y[clade] = leaf_y[clade]
        else:
            child_values = [assign_y(child) for child in clade.clades]
            y[clade] = sum(child_values) / len(child_values)
        return y[clade]

    assign_y(tree.root)
    x = {tree.root: 0.0}

    def assign_x(clade):
        for child in clade.clades:
            x[child] = x[clade] + (child.branch_length if child.branch_length is not None else 1.0)
            assign_x(child)

    assign_x(tree.root)
    names = {leaf.name: y[leaf] for leaf in leaves}
    return x, y, names


def draw_tree(
    ax,
    tree: Tree,
    *,
    label_colors: Mapping[str, str] | None = None,
    stars: set[str] | None = None,
    show_support: bool = True,
    label_size: float = 6.0,
    line_width: float = 0.7,
    label_offset_fraction: float = 0.015,
) -> dict[str, float]:
    label_colors = label_colors or {}
    stars = stars or set()
    x, y, leaf_y = tree_layout(tree)
    max_x = max(x.values()) if x else 1.0

    for clade in tree.find_clades(order="preorder"):
        if clade.clades:
            child_y = [y[child] for child in clade.clades]
            ax.plot([x[clade], x[clade]], [min(child_y), max(child_y)], color="black", lw=line_width)
            if show_support and clade is not tree.root and clade.confidence is not None:
                ax.text(x[clade], y[clade] + 0.12, f"{clade.confidence:g}", fontsize=max(4.5, label_size - 1), ha="center")
        if clade is not tree.root:
            parent = tree.get_path(clade)[-2] if len(tree.get_path(clade)) >= 2 else tree.root
            ax.plot([x[parent], x[clade]], [y[clade], y[clade]], color="black", lw=line_width)

    offset = max(max_x * label_offset_fraction, 0.02)
    for leaf in tree.get_terminals():
        label = leaf.name or ""
        color = label_colors.get(label, "black")
        display = label.replace("_", " ")
        ax.text(max_x + offset, y[leaf], display, fontsize=label_size, color=color, va="center", ha="left", fontstyle="italic")
        if label in stars:
            ax.scatter([max_x - offset * 0.35], [y[leaf]], marker="*", s=(label_size + 1) ** 2, color="#C43C39", linewidths=0, zorder=4)

    ax.set_ylim(-0.8, len(tree.get_terminals()) - 0.2)
    ax.set_xlim(0, max_x + max(offset * 18, max_x * 0.75))
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=max(5, label_size - 1), length=2)
    return leaf_y
