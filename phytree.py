import phytreeviz
from phytreeviz import *

def create_image_of_phytree(filename):
    """
        Creates a phylogenetic tree image from a Newick format file using TreeViz.

        Args:
            filename (str): Path to the Newick file containing the tree structure.

        Returns:
            matplotlib.figure.Figure: A matplotlib figure object containing the rendered tree.
    """

    tv = TreeViz(filename)
    tv.show_branch_length(color="red")
    tv.show_confidence(color="blue")
    tv.show_scale_bar()

    fig = tv.plotfig()
    return fig