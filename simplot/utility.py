
def merged_legend(axes, drawonaxes=None, *args, **kwargs):
    """Merge the legends from the provided list of axes.
    args and kwargs are passed to Axes.legend.
    """
    handles = []
    labels = []
    first = None
    for ax in axes:
        if first is None:
            first = ax
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    if drawonaxes is None:
        drawonaxes = first
    return drawonaxes.legend(handles, labels, *args, **kwargs)