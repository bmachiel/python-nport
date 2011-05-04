""" Defines the LinePlot class.
"""

from __future__ import with_statement

# Standard library imports
import warnings

from numpy import pi, real, imag, ndarray, array, arcsin

from enthought.traits.api import Instance, Property, Float

from enthought.enable.api import ColorTrait, LineStyle

from enthought.chaco.api import Plot, PlotGrid, PlotAxis, AbstractOverlay
from enthought.chaco.api import LinearMapper
from enthought.chaco.api import DataRange1D, ArrayDataSource

from .smith_line_renderer import SmithLineRenderer, SmithCircleRenderer
from enthought.chaco.api import LinearMapper as SmithMapper



class SmithCircle(AbstractOverlay):
    # The color of the axis line.
    line_color = ColorTrait("black")

    # The line thickness (in pixels) of the axis line.
    line_weight = Float(1.0)

    # The dash style of the axis line.
    line_style = LineStyle('solid')
#    def __init__(self, component=None, **kwargs):
#        # Override init so that our component gets set last.  We want the
#        # _component_changed() event handler to get run last.
#        super(self.__class__, self).__init__(**kwargs)
#        if component is not None:
#            self.component = component

    def overlay(self, component, gc, view_bounds=None, mode='normal'):
        """ Draws this component overlaid on another component.

        Overrides AbstractOverlay.
        """
        if not self.visible:
            return
        self._draw_component(gc, view_bounds, mode, component)
        return

#    def _draw_overlay(self, gc, view_bounds=None, mode='normal'):
#        """ Draws the overlay layer of a component.

#        Overrides PlotComponent.
#        """
#        self._draw_component(gc, view_bounds, mode)
#        return

    def _draw_component(self, gc, view_bounds=None, mode='normal', component=None):
        """ Draws the component.

        This method is preserved for backwards compatibility. Overrides
        PlotComponent.
        """
        if not self.visible:
            return

        def center_radius(mapper): #data_low, data_high, screen_low, screen_high):
            map = mapper.map_screen
            return map(0), map(1) - map(0)

        with gc:
            gc.set_stroke_color(self.line_color_)
            gc.set_line_width(self.line_weight)
            gc.set_line_dash(self.line_style_)
            map_x = self.x_mapper.map_screen
            map_y = self.y_mapper.map_screen
            gc.clip_to_rect(map_x(self.x_mapper.range.low),
                            map_y(self.y_mapper.range.low),
                            map_x(self.x_mapper.range.high) -
                              map_x(self.x_mapper.range.low),
                            map_y(self.y_mapper.range.high) -
                              map_y(self.y_mapper.range.low))
                              
            x_center, radius = center_radius(self.x_mapper)
            y_center, radius = center_radius(self.y_mapper)

            # outer ring
            gc.arc(x_center, y_center, radius, 0, 2*pi)
            gc.stroke_path()

            # horizontal axis
            gc.move_to(x_center - radius, y_center)
            gc.line_to(x_center + radius, y_center)
            gc.stroke_path()


class SmithResistanceGrid(AbstractOverlay):
    # The color of the axis line.
    grid_color = ColorTrait("gray")

    # The line thickness (in pixels) of the axis line.
    grid_weight = Float(1.0)

    # The dash style of the axis line.
    grid_style = LineStyle('dot')

    def overlay(self, component, gc, view_bounds=None, mode='normal'):
        """ Draws this component overlaid on another component.

        Overrides AbstractOverlay.
        """
        if not self.visible:
            return
        self._draw_component(gc, view_bounds, mode, component)
        return

    def _draw_component(self, gc, view_bounds=None, mode='normal', component=None):
        """ Draws the component.

        This method is preserved for backwards compatibility. Overrides
        PlotComponent.
        """
        if not self.visible:
            return

        with gc:
            gc.set_stroke_color(self.grid_color_)
            gc.set_line_width(self.grid_weight)
            gc.set_line_dash(self.grid_style_)

            map_x = self.x_mapper.map_screen
            map_y = self.y_mapper.map_screen
            gc.clip_to_rect(map_x(self.x_mapper.range.low),
                            map_y(self.y_mapper.range.low),
                            map_x(self.x_mapper.range.high) -
                              map_x(self.x_mapper.range.low),
                            map_y(self.y_mapper.range.high) -
                              map_y(self.y_mapper.range.low))
                              
            x_center = self.x_mapper.map_screen(0)
            y_center = self.y_mapper.map_screen(0)
            radius = self.y_mapper.map_screen(1) - self.y_mapper.map_screen(0)

            # TODO: adapt to zoom level (fixed number of pixels of spacing)
            for i in range(10):
                r = i * radius / 10
                gc.arc(x_center + r, y_center, radius - r, 0, 2*pi)
                gc.stroke_path()


class SmithReactanceGrid(AbstractOverlay):
    # The color of the axis line.
    grid_color = ColorTrait("gray")

    # The line thickness (in pixels) of the axis line.
    grid_weight = Float(1.0)

    # The dash style of the axis line.
    grid_style = LineStyle('dot')

    def overlay(self, component, gc, view_bounds=None, mode='normal'):
        """ Draws this component overlaid on another component.

        Overrides AbstractOverlay.
        """
        if not self.visible:
            return
        self._draw_component(gc, view_bounds, mode, component)
        return

    def _draw_component(self, gc, view_bounds=None, mode='normal', component=None):
        """ Draws the component.

        This method is preserved for backwards compatibility. Overrides
        PlotComponent.
        """
        if not self.visible:
            return

        with gc:
            gc.set_stroke_color(self.grid_color_)
            gc.set_line_width(self.grid_weight)
            gc.set_line_dash(self.grid_style_)

            map_x = self.x_mapper.map_screen
            map_y = self.y_mapper.map_screen
            gc.clip_to_rect(map_x(self.x_mapper.range.low),
                            map_y(self.y_mapper.range.low),
                            map_x(self.x_mapper.range.high) -
                              map_x(self.x_mapper.range.low),
                            map_y(self.y_mapper.range.high) -
                              map_y(self.y_mapper.range.low))
                              
            x_center = self.x_mapper.map_screen(0)
            y_center = self.y_mapper.map_screen(0)
            rad = self.y_mapper.map_screen(1) - self.y_mapper.map_screen(0)

            # TODO: adapt to zoom level (fixed number of pixels of spacing)
            rad2 = rad**2
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


#-----------------------------------------------------------------------------
# The SmithPlot class
#-----------------------------------------------------------------------------

class SmithPlot(Plot):
    """ Represents a correlated set of data, renderers, and axes in a single
    screen region.

    A Plot can reference an arbitrary amount of data and can have an
    unlimited number of renderers on it, but it has a single X-axis and a
    single Y-axis for all of its associated data. Therefore, there is a single
    range in X and Y, although there can be many different data series. A Plot
    also has a single set of grids and a single background layer for all of its
    renderers.  It cannot be split horizontally or vertically; to do so,
    create a VPlotContainer or HPlotContainer and put the Plots inside those.
    Plots can be overlaid as well; be sure to set the **bgcolor** of the
    overlaying plots to "none" or "transparent".

    A Plot consists of composable sub-plots.  Each of these is created
    or destroyed using the plot() or delplot() methods.  Every time that
    new data is used to drive these sub-plots, it is added to the Plot's
    list of data and data sources.  Data sources are reused whenever
    possible; in order to have the same actual array drive two de-coupled
    data sources, create those data sources before handing them to the Plot.
    """
    
    freq_range = Instance(DataRange1D)

    bgcolor = "transparent"
    border_visible = False
    
    zero_resistance_circle = Instance(SmithCircle)
    resistance_grid = Instance(SmithResistanceGrid)
    reactance_grid = Instance(SmithReactanceGrid)
    
    #------------------------------------------------------------------------
    # Public methods
    #------------------------------------------------------------------------

    def add_smith_plot(self, index_name, value_name, renderer_factory,
                       name=None, origin=None, **kwds):
        raise NotImplementedError

    def add_xy_plot(self, index_name, value_name, renderer_factory, name=None,
        origin=None, **kwds):
        """ Add a BaseXYPlot renderer subclass to this Plot."""
        # doesn't seem to be used
        raise NotImplementedError

    def plot(self, data, name=None, origin=None, **styles):
        """ Adds a new sub-plot using the given data and plot style.

        Parameters
        ==========
        data : string, tuple(string), list(string)
            The data to be plotted. The type of plot and the number of
            arguments determines how the arguments are interpreted:

            one item: (line/scatter)
                The data is treated as the value and self.default_index is
                used as the index.  If **default_index** does not exist, one is
                created from arange(len(*data*))
            two or more items: (line/scatter)
                Interpreted as (index, value1, value2, ...).  Each index,value
                pair forms a new plot of the type specified.
            two items: (cmap_scatter)
                Interpreted as (value, color_values).  Uses **default_index**.
            three or more items: (cmap_scatter)
                Interpreted as (index, val1, color_val1, val2, color_val2, ...)

        name : string
            The name of the plot.  If None, then a default one is created
            (usually "plotNNN").
        origin : string
            Which corner the origin of this plot should occupy:
                "bottom left", "top left", "bottom right", "top right"
        styles : series of keyword arguments
            attributes and values that apply to one or more of the
            plot types requested, e.g.,'line_color' or 'line_width'.

        Examples
        ========
        ::

            plot("my_data", name="myplot", color=lightblue)

            plot(("x-data", "y-data"))

            plot(("x", "y1", "y2", "y3"))

        Returns
        =======
        [renderers] -> list of renderers created in response to this call to plot()
        """
        print "plot"

        if len(data) == 0:
            return

        if isinstance(data, basestring):
            data = (data,)

        # TODO: support lists of plot types
        if name is None:
            name = self._make_new_plot_name()
        if origin is None:
            origin = self.default_origin

        if len(data) == 1:
            # Create the default freq based on the length of the first
            # data series
            value = self._get_or_create_datasource_real(data[0])
            freq = ArrayDataSource(arange(len(value.get_data())),
                                                 sort_order="none")
            #self.freq_range.add(self.default_index)

            if self.default_index is None:
                # Create the default index based the first data series
                value = self._get_or_create_datasource_real(data[0])
                self.default_index = value
                self.index_range.add(self.default_index)
            index = self.default_index
        else:
            index = self._get_or_create_datasource(data[0])
            index_ = self._get_or_create_datasource_real(data[1])
            if self.default_index is None:
                self.default_index = index
            self.index_range.add(index_)
            data = data[1:]

        new_plots = []
        for value_name in data:
            value = self._get_or_create_datasource(value_name)
            value_ = self._get_or_create_datasource_imag(value_name)
            self.value_range.add(value_)
          
            # handle auto-coloring request
            if styles.get("color") == "auto":
                self._auto_color_idx = \
                    (self._auto_color_idx + 1) % len(self.auto_colors)
                styles["color"] = self.auto_colors[self._auto_color_idx]

            imap = SmithMapper(range=self.index_range,
                               stretch_data=self.index_mapper.stretch_data)
            vmap = SmithMapper(range=self.value_range,
                               stretch_data=self.value_mapper.stretch_data)

            plot = SmithLineRenderer(index=index,
                                     value=value,
                                     index_mapper=imap,
                                     value_mapper=vmap,
                                     bgcolor="white",
                                     **styles)

            self.add(plot)
            new_plots.append(plot)

        self.plots[name] = new_plots

        return self.plots[name]

    def plot_circle(self, data, name=None, origin=None, **styles):
        print "plot_circle"
        if len(data) == 0:
            return

        if isinstance(data, basestring):
            data = (data,)

        if name is None:
            name = self._make_new_plot_name()
        if origin is None:
            origin = self.default_origin

        if len(data) == 1:
            if self.default_index is None:
                # Create the default index based on the length of the first
                # data series
                value = self._get_or_create_datasource(data[0])
                self.default_index = ArrayDataSource(arange(len(value.get_data())),
                                                     sort_order="none")
                self.index_range.add(self.default_index)
            index = self.default_index
        else:
            index = self._get_or_create_datasource(data[0])
            index_ = self._get_or_create_datasource_circle_real(data[1])
            if self.default_index is None:
                self.default_index = index
            self.index_range.add(index_)
            data = data[1:]


        new_plots = []
        for value_name in data:
            value = self._get_or_create_datasource(value_name)
            value_ = self._get_or_create_datasource_circle_imag(value_name)
            self.value_range.add(value_)
          
            # handle auto-coloring request
            if styles.get("color") == "auto":
                self._auto_color_idx = \
                    (self._auto_color_idx + 1) % len(self.auto_colors)
                styles["color"] = self.auto_colors[self._auto_color_idx]

            #styles["face_color"] = "green"

            imap = SmithMapper(range=self.index_range,
                               stretch_data=self.index_mapper.stretch_data)
            vmap = SmithMapper(range=self.value_range,
                               stretch_data=self.value_mapper.stretch_data)

            plot = SmithCircleRenderer(index=index,
                                       value=value,
                                       index_mapper=imap,
                                       value_mapper=vmap,
                                       bgcolor="white",
                                       **styles)
                                       
            print "circle:", vmap.range
        
            self.add(plot)
            new_plots.append(plot)

        self.plots[name] = new_plots

        return self.plots[name]

    #------------------------------------------------------------------------
    # Event handlers
    #------------------------------------------------------------------------

    def _bounds_changed(self, old, new):
        print "_bounds_changed"
        super(self.__class__, self)._bounds_changed(old, new)
        if self.value_range.high_setting == 'auto':
            self.index_range.high_setting = 'auto'
        if self.index_range.high_setting == 'auto':
            self.value_range.high_setting = 'auto'

        index_pixels, value_pixels = new[0], new[1]
        index_range = abs(self.index_range.high - self.index_range.low)
        value_range = abs(self.value_range.high - self.value_range.low)
        index_factor = index_pixels / index_range
        value_factor = value_pixels / value_range
        
        if index_factor > value_factor:
            new = self.index_range.low + index_pixels / value_factor
            if new > self.index_range.high:
                self.index_range.high = new
        else:
            new = self.value_range.low + value_pixels / index_factor
            if new > self.value_range.high:
                self.value_range.high = new
        self._update_mappers()


    def _bounds_items_changed(self, event):
        print "_bounds_items_changed"
        super(self.__class__, self)._bounds_items_changed(event)

    def _position_changed(self, old, new):
        print "_position_changed"
        super(self.__class__, self)._position_changed(old, new)

    def _position_items_changed(self, event):
        print "_position_items_changed"
        super(self.__class__, self)._position_items_changed(event)

    def _origin_changed(self):
        print "_origin_changed"
        super(self.__class__, self)._origin_changed()

    def _orientation_changed(self):
        print "_orientation_changed"
        super(self.__class__, self)._orientation_changed()

    def _index_mapper_changed(self, old, new):
        print "_index_mapper_changed"
        super(self.__class__, self)._index_mapper_changed(old, new)

    def _value_mapper_changed(self, old, new):
        print "_value_mapper_changed"
        super(self.__class__, self)._value_mapper_changed(old, new)

    def _x_grid_changed(self, old, new):
        print "_x_grid_changed"
        super(self.__class__, self)._x_grid_changed(old, new)

    def _y_grid_changed(self, old, new):
        print "_y_grid_changed"
        super(self.__class__, self)._y_grid_changed(old, new)

    def _x_axis_changed(self, old, new):
        print "_x_axis_changed"
        super(self.__class__, self)._x_axis_changed(old, new)

    def _y_axis_changed(self, old, new):
        print "_y_axis_changed"
        super(self.__class__, self)._y_axis_changed(old, new)

    def _index_scale_changed(self, old, new):
        print "_index_scale_changed"
        super(self.__class__, self)._index_scale_changed(old, new)

    def _value_scale_changed(self, old, new):
        print "_value_scale_changed"
        super(self.__class__, self)._value_scale_changed(old, new)

    #------------------------------------------------------------------------
    # Private methods
    #------------------------------------------------------------------------

    def _init_components(self):
        # Since this is called after the HasTraits constructor, we have to make
        # sure that we don't blow away any components that the caller may have
        # already set.

        if not self.freq_range:
            self.freq_range = DataRange1D()

        if not self.range2d:
            self.range2d = DataRange2D()

        if not self.index_mapper:
            self.index_mapper = SmithMapper(range=self.range2d.x_range)

        if not self.value_mapper:
            self.value_mapper = SmithMapper(range=self.range2d.y_range)

        # make sure the grid and bgcolor are not the same color

        grid_color = 'lightgray'
#        if color_table[self.bgcolor] == color_table[grid_color]:
#            grid_color = 'white'

        if not self.x_grid:
            self.x_grid = PlotGrid(mapper=self.x_mapper,
                                   orientation="vertical",
                                  line_color=grid_color,
                                  line_style="dot",
                                  component=self)
        if not self.y_grid:
            self.y_grid = PlotGrid(mapper=self.y_mapper,
                                   orientation="horizontal",
                                  line_color=grid_color,
                                  line_style="dot",
                                  component=self)

        if not self.x_axis:
            self.x_axis = PlotAxis(mapper=self.x_mapper,
                                   orientation="bottom",
                                   component=self)

        if not self.y_axis:
            self.y_axis = PlotAxis(mapper=self.y_mapper,
                                   orientation="left",
                                   component=self)

        # Smith axes and grids
        if not self.resistance_grid:
            self.resistance_grid = SmithResistanceGrid(x_mapper=self.x_mapper,
                                                       y_mapper=self.y_mapper,
                                                       component=self)
            self.overlays.append(self.resistance_grid)

        if not self.reactance_grid:
            self.reactance_grid = SmithReactanceGrid(x_mapper=self.x_mapper,
                                                     y_mapper=self.y_mapper,
                                                     component=self)
            self.overlays.append(self.reactance_grid)

        if not self.zero_resistance_circle:
            self.zero_resistance_circle = SmithCircle(x_mapper=self.x_mapper,
                                                      y_mapper=self.y_mapper,
                                                      component=self)
            self.overlays.append(self.zero_resistance_circle)

    def _get_or_create_datasource_circle_real(self, name):
        """ Returns the data source associated with the given name, or creates
        it if it doesn't exist.
        """

        if name + '_real' not in self.datasources:
            data = self.data.get_data(name)

            if type(data) in (list, tuple):
                data = array(data)

            assert isinstance(data, ndarray)

            if isinstance(data, ndarray):
                data = data[:,0].real
                ds = ArrayDataSource(data, sort_order="none")
            elif isinstance(data, AbstractDataSource):
                raise NotImplementedError
                ds = data.real
            else:
                raise ValueError("Couldn't create datasource for data of type " + \
                                 str(type(data)))

            self.datasources[name + '_real'] = ds

        return self.datasources[name + '_real']


    def _get_or_create_datasource_circle_imag(self, name):
        """ Returns the data source associated with the given name, or creates
        it if it doesn't exist.
        """

        if name + '_imag' not in self.datasources:
            data = self.data.get_data(name)

            if type(data) in (list, tuple):
                data = array(data)

            assert isinstance(data, ndarray)

            if isinstance(data, ndarray):
                data = data[:,0].imag
                ds = ArrayDataSource(data, sort_order="none")
            elif isinstance(data, AbstractDataSource):
                raise NotImplementedError
                ds = data.real
            else:
                raise ValueError("Couldn't create datasource for data of type " + \
                                 str(type(data)))

            self.datasources[name + '_imag'] = ds

        return self.datasources[name + '_imag']

    def _get_or_create_datasource_real(self, name):
        """ Returns the data source associated with the given name, or creates
        it if it doesn't exist.
        """

        if name + '_real' not in self.datasources:
            data = self.data.get_data(name)

            if type(data) in (list, tuple):
                data = array(data)

            if isinstance(data, ndarray):
                data = data.real
                assert len(data.shape) == 1
                ds = ArrayDataSource(data, sort_order="none")
            elif isinstance(data, AbstractDataSource):
                raise NotImplementedError
                ds = data.real
            else:
                raise ValueError("Couldn't create datasource for data of type " + \
                                 str(type(data)))

            self.datasources[name + '_real'] = ds

        return self.datasources[name + '_real']

    def _get_or_create_datasource_imag(self, name):
        """ Returns the data source associated with the given name, or creates
        it if it doesn't exist.
        """

        if name + '_imag' not in self.datasources:
            data = self.data.get_data(name)

            if type(data) in (list, tuple):
                data = array(data)

            if isinstance(data, ndarray):
                data = data.imag
                assert len(data.shape) == 1
                ds = ArrayDataSource(data, sort_order="none")
            elif isinstance(data, AbstractDataSource):
                raise NotImplementedError
                ds = data
            else:
                raise ValueError("Couldn't create datasource for data of type " + \
                                 str(type(data)))

            self.datasources[name + '_imag'] = ds

        return self.datasources[name + '_imag']
