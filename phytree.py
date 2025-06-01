import phytreeviz
from phytreeviz import *
import matplotlib.pyplot as plt

def create_image_of_phytree(filename):
    tv = TreeViz(filename)
    tv.show_branch_length(color="red")
    tv.show_confidence(color="blue")
    tv.show_scale_bar()

    fig = tv.plotfig()
    return fig