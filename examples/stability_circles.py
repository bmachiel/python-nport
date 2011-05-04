"""
Draws a static Smith chart plot.
"""

import numpy as np
from nport import touchstone

from enthought.enable.example_support import DemoFrame, demo_main

# Enthought library imports
from enthought.enable.api import Window
from enthought.traits.api import false
from enthought.chaco.api import ArrayPlotData, Plot, OverlayPlotContainer, PlotLabel
from enthought.chaco.tools.api import PanTool, SimpleZoom, DragTool

from smith.smithplot import SmithPlot


class MyFrame(DemoFrame):
    def _create_window(self):
        ts = touchstone.read('../tests/data/deemb_mom.s2p')
        ts = ts.recombine([1,2])

        circles_source = np.vstack(ts.stability_circle_source()).T
        circles_load = np.vstack(ts.stability_circle_load()).T

        container = OverlayPlotContainer(padding=50, fill_padding=True,
                                         bgcolor="lightgray",
                                         use_backbuffer=True)

        self.data = ArrayPlotData(f=ts.freqs,
                                  s11=ts.get_parameter(1 ,1),
                                  s12=ts.get_parameter(1 ,2),
                                  s21=ts.get_parameter(2 ,1),
                                  s22=ts.get_parameter(2 ,2),
                                  circles_source=circles_source,
                                  circles_load=circles_load)
        self.plot = SmithPlot(self.data, title="Smith plot")

        self.plot.plot(("f", "s11"), color="auto", line_width=2.0)
        self.plot.plot(("f", "s22"), color="auto", line_width=2.0)
        #self.plot.plot(("f", "s21"), color="auto", line_width=2.0)
        #self.plot.plot(("f", "s12"), color="auto", line_width=2.0)

        self.plot.plot_circle(("f", "circles_source"), color="auto", line_width=2.0)
        self.plot.plot_circle(("f", "circles_load"), color="auto", line_width=2.0)
       
        container.add(self.plot)

        self.plot.tools.append(PanTool(self.plot))
        zoom = SimpleZoom(self.plot, tool_mode="box", always_on=False)
        self.plot.overlays.append(zoom)

        return Window(self, -1, component=container)


if __name__ == "__main__":
    demo_main(MyFrame, size=(600, 600), title="Simple Smith Chart")
