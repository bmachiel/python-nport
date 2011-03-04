"""
Draws a static Smith chart plot.
"""

from nport import touchstone

from enthought.enable.example_support import DemoFrame, demo_main

# Enthought library imports
from enthought.enable.api import Window
from enthought.traits.api import false
from enthought.chaco.api import OverlayPlotContainer, PlotLabel
from enthought.chaco.tools.api import ZoomTool

# Smith imports
from smith.create_smith_plot import create_smith_plot


class MyFrame(DemoFrame):
    def _create_window(self):
        ts = touchstone.read('tests/data/deemb_mom.s2p')
        freqs = ts.freqs
        s11 = ts.get_parameter(1 ,1)

        container = OverlayPlotContainer(padding=50, fill_padding=True,
                                         bgcolor="lightgray",
                                         use_backbuffer=True)

        plot = create_smith_plot((freqs, s11), color=(0.0, 0.0, 1.0, 1),
                                 width=2.0)
        plot.bgcolor = "white"
        container.add(plot)

        # Add the title at the top
        container.overlays.append(PlotLabel("Smith Chart",
                                  component=container,
                                  font = "swiss 16",
                                  overlay_position="top"))

        return Window(self, -1, component=container)


if __name__ == "__main__":
    demo_main(MyFrame, size=(600, 600), title="Simple Smith Chart")
