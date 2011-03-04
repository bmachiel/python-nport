""" Defines the SmithLineRenderer class.
"""

from __future__ import with_statement

# Major library imports
from numpy import pi, transpose, arcsin

# Enthought library imports
from enthought.traits.api import Float
from enthought.chaco.polar_line_renderer import PolarLineRenderer


class SmithLineRenderer(PolarLineRenderer):
    """ A renderer for Smith plots.
    """
    #------------------------------------------------------------------------
    # Appearance-related traits
    #------------------------------------------------------------------------

    # The width of the origin axis.
    origin_axis_width = 1.0
    # The width of the line.
    line_width = Float(0.5)

    def _draw_background(self, gc, view_bounds=None, mode="default"):
        if self.bgcolor not in ("clear", "transparent", "none"):
            with gc:
                gc.set_antialias(False)
                gc.set_fill_color(self.bgcolor_)
                x_center = self.x + self.width / 2.0
                y_center = self.y + self.height / 2.0

                rad = min(self.width / 2.0, self.height / 2.0)
                gc.arc(x_center, y_center, rad, 0, 2 * pi)
                gc.draw_path()

        # Call the enable _draw_border routine
        if not self.overlay_border and self.border_visible:
            # Tell _draw_border to ignore the self.overlay_border
            self._draw_border(gc, view_bounds, mode, force_draw=True)

        return

    def _draw_default_axes(self, gc):
        if not self.origin_axis_visible:
            return

        with gc:
            gc.set_stroke_color(self.origin_axis_color_)
            gc.set_line_width(self.origin_axis_width)
            gc.set_line_dash(self.line_style_)
            x_center = self.x + self.width / 2.0
            y_center = self.y + self.height / 2.0

            # outer ring
            rad = min(self.width / 2.0, self.height / 2.0)
            gc.arc(x_center, y_center, rad, 0, 2*pi)
            gc.stroke_path()

            # horizontal axis
            gc.move_to(x_center - rad, y_center)
            gc.line_to(x_center + rad, y_center)
            gc.stroke_path()

        return

    def _draw_default_grid(self,gc):
        if not self.grid_visible:
            return

        with gc:
            gc.set_stroke_color(self.origin_axis_color_)
            gc.set_line_width(self.origin_axis_width)
            gc.set_line_dash(self.grid_style_)
            x_data,y_data = transpose(self._cached_data_pts)
            x_center = self.x + self.width / 2.0
            y_center = self.y + self.height / 2.0
            rad = min(self.width / 2.0, self.height / 2.0)
            rad2 = rad**2
            for i in range(10):
                r = i * rad / 10
                gc.arc(x_center + r, y_center, rad - r, 0, 2*pi)
                gc.stroke_path()

            for i in range(6):
                r1 = i * rad / 5
                tmp = (rad2 - r1**2) / (rad2 + r1**2)
                gc.arc(x_center + rad, y_center + r1, r1, pi - arcsin(tmp), 1.5*pi)
                gc.stroke_path()
                gc.arc(x_center + rad, y_center - r1, r1, 0.5*pi, pi + arcsin(tmp))
                gc.stroke_path()

            stoprad = 2 * rad / 10
            stoprad2 = stoprad**2
            for i in range(6):
                r2 = 7 * rad / (i + 1)
                tmp = (rad2 - r2**2) / (rad2 + r2**2)
                tmp2 = (stoprad2 - r2**2) / (stoprad2 + r2**2)
                gc.arc(x_center + rad, y_center + r2, r2,
                       pi - arcsin(tmp), pi - arcsin(tmp2))
                gc.stroke_path()
                gc.arc(x_center + rad, y_center - r2, r2,
                       pi + arcsin(tmp2), pi + arcsin(tmp))
                gc.stroke_path()

        return

