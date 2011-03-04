"""
Contains a convenience function to create a ready-made Smith plot
"""

# Major library imports
from numpy import array, ndarray, transpose, real, imag

from enthought.chaco.array_data_source import ArrayDataSource
from enthought.chaco.data_range_1d import DataRange1D

# Local relative imports
from .smith_mapper import SmithMapper
from .smith_line_renderer import SmithLineRenderer


def create_smith_plot(data, orientation='h', color='black', width=1.0,
                      dash="solid", grid="dot", value_mapper_class=SmithMapper):
    if (type(data) != ndarray) and (len(data) == 2):
        data = transpose(array(data))

    freqs, values = transpose(data)
    index_data = real(values)
    value_data = imag(values)

    index = ArrayDataSource(index_data, sort_order='ascending')
    # Typically the value data is unsorted
    value = ArrayDataSource(value_data)

    index_range = DataRange1D()
    index_range.add(index)
    index_mapper = SmithMapper(range=index_range)

    value_range = DataRange1D()
    value_range.add(value)
    value_mapper = value_mapper_class(range=value_range)

    plot = SmithLineRenderer(index=index, value=value,
                             index_mapper = index_mapper,
                             value_mapper = value_mapper,
                             orientation = orientation,
                             color = color,
                             line_width = width,
                             line_style = dash,
                             grid_style = grid)
    return plot
