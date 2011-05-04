""" Defines the SmithLineRenderer class.
"""

from __future__ import with_statement

# Major library imports
import numpy as np
from numpy import array, pi, transpose, arcsin, real, imag, hstack

# Enthought library imports
from enthought.enable.api import LineStyle, black_color_trait, \
                                  transparent_color_trait
from enthought.kiva.agg import points_in_polygon
from enthought.traits.api import Enum, Float
from enthought.chaco.api import LinePlot, BaseXYPlot


class SmithLineRenderer(LinePlot):
#    def _gather_points(self):
#        print "SmithLineRenderer._gather_points"
##        self._cache_valid = False
#        print "index", self.index.get_data()
#        print "value", self.value.get_data()
#        super(self.__class__, self)._gather_points()

    def _render(self, gc, points):
        print "SmithLineRenderer._render"
        super(self.__class__, self)._render(gc, points)

    def _gather_points(self):
        """ Collects the data points that are within the bounds of the plot and
        caches them.
        """
        if self._cache_valid:
            return

        data = self.value.get_data()
        x = data.real.astype(float)
        y = data.imag.astype(float)
        
        if not self.index or not self.value:
            return

        if len(x) == 0:
            self._cached_data_pts = []
            self._cache_valid = True
            return

        points = [np.transpose(np.array((x, y)))]
        self._cached_data_pts = points
        self._cache_valid = True

#    def _render(self, gc, circles):
#        """ Renders a circle."""
#        map_x = self.index_mapper.map_screen
#        map_y = self.value_mapper.map_screen

#        for circle in circles:
#            x_center = map_x(circle[0])
#            y_center = map_y(circle[1])
#            radius = map_x(circle[2]) - map_x(0)
#            with gc:
#                gc.clip_to_rect(map_x(self.x_mapper.range.low),
#                                map_y(self.y_mapper.range.low),
#                                map_x(self.x_mapper.range.high) -
#                                  map_x(self.x_mapper.range.low),
#                                map_y(self.y_mapper.range.high) -
#                                  map_y(self.y_mapper.range.low))

#                gc.set_stroke_color(self.color_)
#                gc.set_line_width(self.line_width)
#                gc.set_line_dash(self.line_style_)
#                gc.set_fill_color(self.face_color_)

#                gc.arc(x_center, y_center, radius, 0, 2*pi)
#                gc.draw_path()


class SmithCircleRenderer(BaseXYPlot):
    """ Plots a circle in dataspace.

    If you don't want the edge of the polygon to be drawn, set **edge_color**
    to transparent; don't try to do this by setting **edge_width** to 0. In
    some drawing systems, such as PostScript, a line width of 0 means to make
    the line as small as possible while still putting ink on the page.
    """

    # The color of the line on the edge of the polygon.
    color = black_color_trait

    # The thickness of the edge of the polygon.
    line_width = Float(1.0)

    # The line dash style for the edge of the polygon.
    line_style = LineStyle

    # The color of the face of the polygon.
    face_color = transparent_color_trait

    # Override the hittest_type trait inherited from BaseXYPlot
    hittest_type = Enum("poly", "point", "line")

    #### Private 'BaseXYPlot' interface ########################################

#    def __init__(self, *args, **kw):
#        super(self.__class__, self).__init__(*args, **kw)

    def _downsample(self):
        #return self.map_screen(self._cached_data_pts)
        return self._cached_data_pts

    def _draw_plot(self, *args, **kw):
        """ Draws the 'plot' layer.
        """
        # Simple compatibility with new-style rendering loop
        return self._draw_component(*args, **kw)

    def _draw_component(self, gc, view_bounds=None, mode='normal'):
        """ Renders the component. 
        """
        self._gather_points()
        self._render(gc, self._cached_data_pts)

    def _gather_points(self):
        """ Collects the data points that are within the bounds of the plot and
        caches them.
        """
        if self._cache_valid:
            return

        data = transpose(self.value.get_data())
        x_center = data[0].real.astype(float)
        y_center = data[0].imag.astype(float)
        radius = data[1].astype(float)
        
        if not self.index or not self.value:
            return

        if len(x_center) == 0 or len(radius) == 0 or len(x_center) != len(radius):
            self._cached_data_pts = []
            self._cache_valid = True
            return

        circles = np.transpose(np.array((x_center, y_center, radius)))
        self._cached_data_pts = circles
        self._cache_valid = True

    def _render(self, gc, circles):
        """ Renders a circle."""
        map_x = self.index_mapper.map_screen
        map_y = self.value_mapper.map_screen

        for circle in circles:
            x_center = map_x(circle[0])
            y_center = map_y(circle[1])
            radius = map_x(circle[2]) - map_x(0)
            with gc:
                gc.clip_to_rect(map_x(self.x_mapper.range.low),
                                map_y(self.y_mapper.range.low),
                                map_x(self.x_mapper.range.high) -
                                  map_x(self.x_mapper.range.low),
                                map_y(self.y_mapper.range.high) -
                                  map_y(self.y_mapper.range.low))

                gc.set_stroke_color(self.color_)
                gc.set_line_width(self.line_width)
                gc.set_line_dash(self.line_style_)
                gc.set_fill_color(self.face_color_)

                gc.arc(x_center, y_center, radius, 0, 2*pi)
                gc.draw_path()

    def _render_icon(self, gc, x, y, width, height):
        """ Renders a representation of this plot as an icon into the box
        defined by the parameters.

        Used by the legend.
        """
        with gc:
            gc.set_stroke_color(self.edge_color_)
            gc.set_line_width(self.line_width)
            gc.set_fill_color(self.face_color_)
            if hasattr(self, 'line_style_'):
                gc.set_line_dash(self.line_style_)
            gc.draw_rect((x,y,width,height))
        return

    #------------------------------------------------------------------------
    # AbstractPlotRenderer interface
    #------------------------------------------------------------------------

#    def map_screen(self, data_array):
#        """ Maps an array of data points into screen space and returns it as
#        an array. 
        
#        Implements the AbstractPlotRenderer interface.
#        """
#        print "map_screen"
#        # data_array is Nx2 array
#        if len(data_array) == 0:
#            return []
#        x_ary, y_ary, r_ary = transpose(data_array)
#        sx = self.index_mapper.map_screen(x_ary)
#        sy = self.value_mapper.map_screen(y_ary)
#        sr = (self.index_mapper.map_screen(r_ary) -
#              self.index_mapper.map_screen(0))
#        for i in range(len(sx)):
#            print "center = ({0}, {1}), radius = {2}".format(x_ary[i], y_ary[i],
#                                                             r_ary[i])
#        if self.orientation == "h":
#            return transpose(array((sx, sy, sr)))
#        else:
#            return transpose(array((sy, sx, sr)))

    def hittest(self, screen_pt, threshold=7.0, return_distance=False):
        """ Performs point-in-polygon testing or point/line proximity testing.
        If self.hittest_type is "line" or "point", then behaves like the
        parent class BaseXYPlot.hittest().

        If self.hittest_type is "poly", then returns True if the given
        point is inside the polygon, and False otherwise.
        """
        if self.hittest_type in ("line", "point"):
            return BaseXYPlot.hittest(self, screen_pt, threshold, return_distance)

        data_pt = self.map_data(screen_pt, all_values=True)
        index = self.index.get_data()
        value = self.value.get_data()
        poly = np.vstack((index,value)).T
        if points_in_polygon([data_pt], poly)[0] == 1:
            return True
        else:
            return False

    #------------------------------------------------------------------------
    # Event handlers
    #------------------------------------------------------------------------

    def _alpha_changed(self):
        self.color_ = self.color_[0:3] + (self.alpha,)
        self.face_color_ = self.face_color_[0:3] + (self.alpha,)
        self.invalidate_draw()
        self.request_redraw()