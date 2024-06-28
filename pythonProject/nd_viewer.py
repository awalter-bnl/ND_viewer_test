import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import xarray as xr


class LinkableSlider(mpl.widgets.Slider):
    """A child of mpl.widgets.Slider which allows linking sliders to each other.

    Essentially this adds a new attribute 'self.linked_sliders' which is an
    updatable list of sliders that should change with this slider. It also adds
    attributes that lists the 'plots', cursors and plot_parameters associated
    with this slider as well. The update function that is attached to the slider
     should take care of updating all linked sliders as well.

    Attributes:
    -----------
    linked_sliders: [LinkableSlider, ....].
        A list of LinkableSlider or matplotlib.widgets.Slider objects that
        should be updated, and have the associated update functions called, when
        this LinkableSlider is updated.
    See mpl.widgets.Slider documentation for additional attributes.

    Methods:
    --------
    set_val(val):
        In addition to calling the parent.set_val method it also ensures that
        each of the LinkableSliders/  mpl.widgets.Slider objects in the
        linked_sliders attribute list.
    """

    def __init__(self, *args, plots, cursors, plot_parameters, **kwargs):
        """The init method.

        Parameters:
        -----------
        plots: dict
            Dictionary mapping tuples of plotted dimension names to
            matplotlib.AxesImage or matplotlib.line2D objects for each of the
            2D(1D) slice plots associated with this slider.
        cursors: dict
            Dictionary mapping tuples of plotted dimension names to
            matplotlib.line2D objects for each of the 1D slice plots associated
            with this slider.
        plot_parameters: dict
            Dictionary mapping tuples of plotted dimension names to a dictionary
             of the form :
                {'transpose': Bool,'flip_dims':List} where transpose indicates
                 if an axes should be transposed and 'flip_dims' is a list that
                 indicates which axes need to be flipped.
        *args: list
            A list of arguments passed to the parent init class
        **kwargs:

        Passes *args and **kwargs to the parent init class and adds the new
        attribute 'linked_sliders'.
        """
        self.linked_sliders = []
        self.plots = plots
        self.cursors = cursors
        self.plot_parameters = plot_parameters
        super().__init__(*args, **kwargs)


def multislice(*view_args, controls={}, buttons={}, fig_scale=1,
               close_figures=True, **kwargs):
    '''Allows to plot multiple linked 'multislice views' in the same figure.

    This function allows for multiple, linked, 'multislice views', generated
    using 'multislice_view' to be plotted in the same figure. It creates a
    global figure, with subfigures for each 'view', and then creates each of the
     views keeping track of the 'controls' object returned so that they are all
     linked. For a description of the plots displayed please see the description
      of 'multiuslice_view'.

    Parameters:
    -----------
    *view_args :
        Has the structure:

        .. code-block:: python
            data1, (main_dims1, title1,)
            data2, (main_dims2, title2,)
            ...,
            dataN, (main_dimsN, titleN)

        where data* is an xarray that requires a multislice view, main_dims* is
        an optional list of dimensions to use for axes in the prior data*
        multislice view and title* is the optional title to use for the prior
        data* multislice view.
    controls : {str:[LinkedSlider,...](,str:[LinkedSlider,...],....)}.
        A dictionary mapping dimension names to the sliders that update those
        dimensions.
    buttons : {str:{'left':[mpl.widgets.Button,...],
                    'left':[mpl.widgets.Button,...]},....}.
        A dictionary mapping dimension names to the buttons that update those
        dimensions.
    fig_scale : float, optional.
        An optional float that will scale the size of the figure including any
        text size set using fontsize, ignored if 'fig' is passed in as that
        determines the figure size. Default size is made so that the full height
         of the figure can be displayed in jupyter with the browser scale set to
          100% on a standard 14" laptop monitor.
    **kwargs :
        kwargs passed to each call of 'multislice_view', see it's description
        for details.
    close_figures: Bool, optional.
        A boolean to indicate if the plot should call matplotlib.pyplot.close()
        before continuing. This has proven useful to avoid a 'too many figures
        open' warning in Jupyter notebooks.

    Returns:
    --------
    fig: matplotlib.pyplot.figure.
        A matplotlib.pyplot.figure object that the plot was generated on, is
        passed through if 'fig' is provided and otherwise it is the created
        figure object.
    controls : {str:[LinkedSlider,...](,str:[LinkedSlider,...],....)}.
        A dictionary mapping dimension names to the sliders that update those
        dimensions, it is an updated version of the object that is passed in
        with any additional sliders required for this view added.
    buttons : {str:{'left':[mpl.widgets.Button,...],
                    'left':[mpl.widgets.Button,...]},....}.
        A dictionary mapping dimension names to the buttons that update those
        dimensions, it is an updated version of the object that is passed in
        with any additional buttons required for this view added.

    '''

    def _parse_args(view_args):
        '''Parses the arguments passed in to a list of 3 value lists associated
        with each xarray listed in args.

        This will loop through all the items in *args and create a list of 3
        element lists, one for each xarray in *args with the other two elements
        giving the list of main_dims and the 3rd giving the title for the xarray
        multislice view.

        Parameters
        ----------
        *view_args:
            Has the structure:

            .. code-block:: python
                data1, (main_dims1, title1,)
                data2, (main_dims2, title2,)
                ...,
                dataN, (main_dimsN, titleN)

            where data* is an xarray that requires a multislice view, main_dims*
            is an optional list of dimensions to use for axes in the prior
            data* multislice view and title* is the optional title to use for
            the prior data* multislice view.

        Returns
        -------
        split: [[data1, main_dims1, title1], ... , [dataN, main_dimsN, titleN]]
            A list of 3 element lists for each xarray for which a multislice
            view is required.

        '''

        split = []
        temp = []
        position = 0
        for i, item in enumerate(view_args):
            if type(item) is xr.DataArray:  # if the item is a data array
                if temp: split.append(temp)
                position = 0
                temp = [item, None, f'Multislice View {i + 1}']
            elif type(item) is list:
                temp[1] = item
            elif type(item) is str:
                temp[2] = item
            else:
                raise ValueError(f'Expected items in args to contain an '
                                 f'xarray.DataArray, a list or a str instead '
                                 f'got {item} in {view_args}')

        if temp: split.append(temp)  # add the last group

        return split

    # close figures to prevent to many being open if asked.
    if close_figures:
        plt.close()

    view_data = _parse_args(view_args)

    # Determine Figure size and create subfigures for each 'view'.
    fig_delta = 0.12
    num_controls = max([max(3, len(data.dims) - 3)
                        for (data, main_dims, title) in view_data])
    fig = plt.figure(figsize=[len(view_data) * 40 * fig_delta * fig_scale,
                              (20 * fig_delta + 22 * fig_delta +
                               num_controls * fig_delta) * fig_scale])
    subfigs = fig.subfigures(1, len(view_data))
    if len(view_data) == 1:
        subfigs = np.array([subfigs])
    # Required workaround to an mpl bug where the default facecolor of
    # subfigures is not transparent
    for subfig in subfigs.flat:
        subfig.set_facecolor((1, 1, 1, 0.1))

        # generate each sub-view
    for i, (data, main_dims, title) in enumerate(view_data):
        temp_fig, controls, buttons = multislice_view(data, main_dims=main_dims,
                                                      fig=subfigs[i],
                                                      controls=controls,
                                                      buttons=buttons,
                                                      figure_title=title,
                                                      close_figures=False,
                                                      **kwargs)

    return fig, controls, buttons


def multislice_view(data, main_dims=None, fig=None, controls={}, buttons={},
                    figure_title='Multislice View', fontsize=8,
                    cursor_colours=['lightgrey', 'orange', 'yellow'],
                    close_figures=True, fig_scale=1):
    '''Displays 2D and 1D 'slices' through 'N-D' datasets.

    This function is designed to display 1(3) 2D 'slices' and the corresponding
    1D 'slices' from an 'ND' dataset. It displays a single 2D/two 1D slices if
    the 'main_dims' parameter has 2 values and three 2D/three 1D slices if
    'main_dims' has 3 values. It will automatically detect all axes in the
    'data' xarray and provide controls for adjusting which value to use for
    these when performing the slicing. Note that the higher level 'multislice'
    function can be used if you have multiple datasets that need to be
    controlled simultaneously.

    This function uses the packages matplotlib, mpl_interactions for the
    plotting and matplotlib.widgets objects for the controls and so should be
    compatible with any plotting back-end that they are.

    Parameters
    ----------
    data: xarray.DataArray.
        An N-dimensional xarray.DataArray object that contains the data to be
        plotted.
    main_dims: [str,str(,str)], optional.
        A two(three) element list or tuple that indicates which dimensions from
        data should be used as the axes for the 2D slice(s). The first two
        dimensions are used for the axes of the 'central' 2D slice while the
        optional third dimension will be combined with the first two to display
        two more slices. Each of these dimensions will have their own associated
         1D slice. If they are not provided it takes the first two dimensions of
          data for the 'central' axis, and if a third dimension exists this will
           be used for the additional two axes.
    fig: matplotlib.figure.Figure, optional.
        A matplotlib.figure.Figure object to be used to generate the plot on, if
         not given it generates it's own. This is provided to allow this
         function to be used to generate a sub-figure in a larger plot.
    controls : {str:[LinkedSlider,...](,str:[LinkedSlider,...],....)}.
        A dictionary mapping dimension names to the sliders that update those
        dimensions.
    buttons : {str:{'left':[mpl.widgets.Button,...],
                    'left':[mpl.widgets.Button,...]},....}.
        A dictionary mapping dimension names to the buttons that update those
        dimensions.
    figure_title: str, optional
        An optional figure title string, if not provided it will use a generic
        title.
    fontsize: float or string, optional.
        Tick label font size in points or as a string (e.g., 'large') passed to
        the matplotlib labels and ticks. Note
        that the text is scaled if the optional fig_scale is not 1.
    cursor_colours: [matplotlib color format, matplotlib color format,
                     matplotlib color format], optional.
        A three element list/tuple indicating the colour of the cursors for each
         axis value, if plotting 1 2D slice then the third element is ignored.
    close_figures: Bool, optional.
        A boolean to indicate if the plot should call matplotlib.pyplot.close()
        before continuing. This has proven useful to avoid a 'too many figures
        open' warning in Jupyter notebooks.
    fig_scale: float, optional.
        An optional float that will scale the size of the figure including any
        text size set using fontsize, ignored if 'fig' is passed in as that
        determines the figure size. Default size is made so that the full height
         of the figure can be displayed in jupyter with the browser scale set to
          100% on a standard 14" laptop monitor.

    Returns
    -------
    fig: matplotlib.pyplot.figure.
        A matplotlib.pyplot.figure object that the plot was generated on, is
        passed through if 'fig' is provided and otherwise it is the created
        figure object.
    controls : {str:[LinkedSlider,...](,str:[LinkedSlider,...],....)}.
        A dictionary mapping dimension names to the sliders that update those
        dimensions, it is an updated version of the object that is passed in
        with any additional sliders required for this view added.
    buttons : {str:{'left':[mpl.widgets.Button,...],
                    'left':[mpl.widgets.Button,...]},....}.
        A dictionary mapping dimension names to the buttons that update those
        dimensions, it is an updated version of the object that is passed in
        with any additional buttons required for this view added.

    '''

    def _generate_layout(main_dims, extra_dims, fig=None):
        '''Generates the figure layouts and axes.

        sub-function that generates the layout of the plot and returns
        dictionaries mapping dimension names to matplotlib.axes for the 2D
        plots, 1D plots and the control widgets.

        main_dims: [str,str(,str)].
            A two(three) element list or tuple that indicates which dimensions
            from data should be used as the axes for the 2D slice(s). The first
            two dimensions are used for the axes of the 'central' 2D slice while
            the optional third dimension will be combined with the first two to
            display two more slices. Each of these dimensions will have their
            own associated 1D slice.
        extra_dims : [str, str, ...].
            A list of strings containing all of the dimensions, not included in
            main_dims, that the layout should provide controls axes for.
        fig: matplotlib.figure.Figure, optional.
            A matplotlib.figure.Figure object to be used to generate the plot
            on, if not given it generates it's own. This is provided to allow
            this function to be used to generate a sub-figure in a larger plot.

        Returns
        -------
        fig: matplotlib.figure.Figure.
            The matplotlib.figure.Figure object that the layout is made on.
        axes2D,axes1D: Dicts.
            Dictionary mapping tuples of plotted dimension names to
            matplotlib.axes objects for each of the 2D(1D) slice plots.
        axesC: {str:{'slider':matplotlib.axes,'left button':matplotlib.axes,
                     'right button':matplotlib.axes}}.
            Dictionary mapping dimension names to a dictionary mapping the keys
            'slider', 'left button', 'right button' to matplotlib.axes objects
            to be used for the widget controlling that dimension.
        '''
        if len(main_dims) == 3:
            indices_map2D = {(main_dims[2], main_dims[1]): [0, 0],
                             (main_dims[0], main_dims[2]): [1, 1],
                             (main_dims[0], main_dims[1]): [0, 1]}
            layout_map = {
                (main_dims[2], main_dims[1]):
                    {'left': 0.19, 'right': 0.99, 'bottom': 0.01, 'top': 0.76},
                (main_dims[0], main_dims[2]):
                    {'left': 0.01, 'right': 0.81, 'bottom': 0.24, 'top': 0.99},
                (main_dims[0], main_dims[1]):
                    {'left': 0.01, 'right': 0.81, 'bottom': 0.01, 'top': 0.76}}
            indices1D = [1, 0]  # gives the 1D plot grid indices.
            control_index = 2  # gives the row index for the controls

        elif len(main_dims) == 2:
            indices_map2D = {(main_dims[0], main_dims[1]): [0, 1]}
            layout_map = {
                (main_dims[0], main_dims[1]):
                    {'left': 0.01, 'right': 0.81, 'bottom': 0.01, 'top': 0.76}}
            indices1D = [0, 0]  # gives the 1D plot grid indices.
            control_index = 1  # gives the row index for the controls

        num_controls = max(len(main_dims), len(extra_dims))
        if not fig:
            fig_delta = 0.12
            fig = plt.figure(figsize=[40 * fig_delta * fig_scale,
                                      (20 * fig_delta + 22 * fig_delta +
                                       num_controls * fig_delta) * fig_scale])

        fig.suptitle(figure_title)
        height_ratios = [len(main_dims) - 1, 0.05 * num_controls,
                         abs(len(main_dims) - 3)]
        subfigs = fig.subfigures(3, 1, height_ratios=height_ratios)
        # Required workaround to an mpl bug where the default facecolor of
        # subfigures is not transparent
        for subfig in subfigs.flat:
            subfig.set_facecolor((1, 1, 1, 0.1))

        plot_subfigs = subfigs[0].subfigures(
            len(main_dims) - 1, 2, width_ratios=[1, 1],
            squeeze=False)  # splits plot subfigs
        # Required workaround to an mpl bug where the default facecolor of
        # subfigures is not transparent
        for subfig in plot_subfigs.flat:
            subfig.set_facecolor((1, 1, 1, 0.1))

            # 2D layout
        axes2D = {}
        for axes_dims, indices2D in indices_map2D.items():
            axes2D[axes_dims] = plot_subfigs[
                indices2D[0],
                indices2D[1]].subplots(1, 1, gridspec_kw=layout_map[axes_dims])
            if len(main_dims) == 3:
                if axes_dims == (main_dims[2], main_dims[1]):
                    axes2D[axes_dims].set_xlabel(data.dims[2],
                                                 fontsize=fontsize * fig_scale)
                    axes2D[axes_dims].xaxis.set_label_position('top')
                    axes2D[axes_dims].set_ylabel(data.dims[1],
                                                 fontsize=fontsize * fig_scale)
                    axes2D[axes_dims].yaxis.set_label_position('left')
                    axes2D[axes_dims].invert_xaxis()
                    axes2D[axes_dims].tick_params(top=True, labeltop=True,
                                                  bottom=True, labelbottom=True,
                                                  left=True, labelleft=True,
                                                  right=True, labelright=False,
                                                  labelsize=fontsize*fig_scale)

                elif axes_dims == (main_dims[0], main_dims[2]):
                    axes2D[axes_dims].set_xlabel(data.dims[0],
                                                 fontsize=fontsize * fig_scale)
                    axes2D[axes_dims].set_ylabel(data.dims[2],
                                                 fontsize=fontsize * fig_scale)
                    axes2D[axes_dims].yaxis.set_label_position('right')
                    axes2D[axes_dims].tick_params(top=True, labeltop=False,
                                                  bottom=True, labelbottom=True,
                                                  left=True, labelleft=False,
                                                  right=True, labelright=True,
                                                  labelsize=fontsize*fig_scale)
                elif axes_dims == (main_dims[0], main_dims[1]):
                    axes2D[axes_dims].set_xlabel(data.dims[0],
                                                 fontsize=fontsize * fig_scale)
                    axes2D[axes_dims].xaxis.set_label_position('top')
                    axes2D[axes_dims].set_ylabel(data.dims[1],
                                                 fontsize=fontsize * fig_scale)
                    axes2D[axes_dims].yaxis.set_label_position('right')
                    axes2D[axes_dims].tick_params(top=True, labeltop=True,
                                                  bottom=True, labelbottom=True,
                                                  left=True, labelleft=False,
                                                  right=True, labelright=True,
                                                  labelsize=fontsize*fig_scale)
            else:
                if axes_dims == (main_dims[0], main_dims[1]):
                    axes2D[axes_dims].set_xlabel(data.dims[0],
                                                 fontsize=fontsize * fig_scale)
                    axes2D[axes_dims].xaxis.set_label_position('top')
                    axes2D[axes_dims].set_ylabel(data.dims[1],
                                                 fontsize=fontsize * fig_scale)
                    axes2D[axes_dims].yaxis.set_label_position('right')
                    axes2D[axes_dims].tick_params(top=True, labeltop=True,
                                                  bottom=True, labelbottom=False,
                                                  left=True, labelleft=False,
                                                  right=True, labelright=True,
                                                  labelsize=fontsize*fig_scale)

        # 1D layout
        if len(main_dims) == 3:
            layout = {'left': 0.19, 'right': 0.99, 'bottom': 0.12, 'top': 0.99,
                      'hspace': 0.65}
        else:
            layout = {'left': 0.19, 'right': 0.99, 'bottom': 0.15, 'top': 0.76,
                      'hspace': 0.65}
        axes1D = {}
        for i, subplot in enumerate(plot_subfigs[
            indices1D[0],
            indices1D[1]].subplots(len(main_dims), 1, gridspec_kw=layout)):
            axes1D[(main_dims[i],)] = subplot
            axes1D[(main_dims[i],)].set_xlabel(data.dims[i], labelpad=1.0,
                                               fontsize=fontsize * fig_scale)
            axes1D[(main_dims[i],)].set_ylabel('Arb. units', labelpad=1.0,
                                               fontsize=fontsize * fig_scale)
            axes1D[(main_dims[i],)].tick_params(top=False, labeltop=False,
                                                bottom=True, labelbottom=True,
                                                left=True, labelleft=True,
                                                right=False, labelright=False,
                                                labelsize=fontsize * fig_scale)

        # Control layout
        axesC = {}
        controls_subfigs = subfigs[1].subfigures(1, 2)
        # Required workaround to an mpl bug where the default facecolor of
        # subfigures is not transparent
        for subfig in controls_subfigs.flat:
            subfig.set_facecolor((1, 1, 1, 0.1))

        main_subfigs = controls_subfigs[1].subfigures(num_controls, 4,
                                                      width_ratios=[14, 1, 1,
                                                                    2])
        for i, dim in enumerate(main_dims):
            axesC[dim] = {'slider': main_subfigs[i][0].add_axes([0.1, 0, 0.65,
                                                                 1]),
                          'left_button': main_subfigs[i][1].add_axes([0, 0, 1,
                                                                      1]),
                          'right_button': main_subfigs[i][2].add_axes([0, 0, 1,
                                                                       1])}

        extra_subfigs = controls_subfigs[0].subfigures(num_controls, 4,
                                                       width_ratios=[2, 14, 1,
                                                                     1])
        for i, dim in enumerate(extra_dims):
            axesC[dim] = {'slider': extra_subfigs[i][1].add_axes([0.1, 0, 0.65,
                                                                  1]),
                          'left_button': extra_subfigs[i][2].add_axes([0, 0, 1,
                                                                       1]),
                          'right_button': extra_subfigs[i][3].add_axes([0, 0, 1,
                                                                        1])}

        return fig, axes2D, axes1D, axesC

    def _slider_update_factory(slider, dim):
        '''A function generator for the slider widgets.

        XXXXX

        Parameters
        ----------
        slider: mpl.widgets.Slider,
            The mpl.widgets.Slider that the update function is meant to work on.
        dim : str.
            The name of the dimension that this slider adjusts

        Returns
        -------
        _update : function,
            The update function that can be used by the slider widget.
        '''

        def _update(val):
            '''An specific update function to be used by widgets.

            Parameters
            ----------
            val : float.
                The value that the slider returns.
            '''

            def _slider_update(_slider, val):
                # Update the associated plots
                for axis_dims, plot in _slider.plots.items():
                    if dim not in axis_dims:
                        indexers = {d: sliders[0].val
                                    for d, sliders in controls.items()
                                    if d not in axis_dims}
                        _slice = data.sel(indexers, method="nearest")

                        if len(axis_dims) == 2:  # 2Dplot
                            if _slider.plot_parameters[axis_dims]['transpose']:
                                _slice = _slice.T
                            for flip_dim in _slider.plot_parameters[axis_dims][
                                'flip_dims']:
                                _slice = _slice.isel({flip_dim: slice(None,
                                                                      None,
                                                                      -1)})
                            plot.set_data(_slice)
                            plot.norm.autoscale([float(_slice.min()),
                                                 float(_slice.max())])
                            plot.axes.draw_artist(plot)
                        elif len(axis_dims) == 1:  # 1Dplot
                            for flip_dim in _slider.plot_parameters[axis_dims][
                                'flip_dims']:
                                _slice = _slice.isel({flip_dim: slice(None,
                                                                      None,
                                                                      -1)})
                            plot.set_ydata(_slice)
                            plot.axes.set_ylim(float(_slice.min()),
                                               float(_slice.max()))
                            plot.axes.draw_artist(plot)
                        else:
                            raise ValueError(
                                f'In a multislice_view._slider_update_factory.'
                                f'_update call Expected the keys in '
                                f'_slider.plots to be a 1D or 2D list of '
                                f'dimension names but instead got {axis_dims}.')

                # Update the associated Cursors
                for axis_dims, axis_cursors in _slider.cursors.items():
                    if axis_dims[0] == dim:
                        for x_cursor in axis_cursors['x']:
                            x_cursor.set_xdata([val, val])
                            x_cursor.axes.draw_artist(x_cursor)
                    if len(axis_dims) == 2 and axis_dims[1] == dim:
                        for y_cursor in axis_cursors['y']:
                            y_cursor.set_ydata([val, val])
                            y_cursor.axes.draw_artist(y_cursor)
            _slider_update(slider, val)

            # update the linked sliders
            for linked_slider in slider.linked_sliders:
                _slider_update(linked_slider, val)
                eventson_state = linked_slider.eventson
                linked_slider.eventson = False
                linked_slider.set_val(val)
                linked_slider.eventson = eventson_state

        return _update

    def _button_update_factory(slider, direction='forward'):
        '''Returns an update function for forward/back buttons/

        This is a wrapper function that will return an update function that
        increases/decreases the value of 'slider' by one on clicking.

        Parameters
        ----------
        slider: mpl.widgets.Slider,
            The mpl.widgets.Slider object that should be updated on button
            click.
        direction: str, optional.
            The 'direction' that the button should move along the dimensions of
            'slider', it should be 'forward' or 'backward'.

        Returns
        -------
        _update : function,
            The update function that can be used by the button widget.
        '''
        if direction not in ['forward', 'backward']:
            raise ValueError('In multislice_layout._update_button_wrapper the '
                             '"direction" parameter is not "forward" or '
                             '"backward" as required. direction = {direction}')

        def _update(*args, **kwargs):
            index_step = 1 if direction == 'forward' else -1
            current_index = np.abs(slider.valstep.values - slider.val).argmin()
            slider.set_val(slider.valstep.values[current_index + index_step])

        return _update

    if not controls:  # Not sure why I need to do this but the variable controls
        # acts like a global if passed in
        controls = {}

    ### FUNCTION STARTS HERE
    # Add some ValueError checks here.
    if type(data) not in [xr.core.dataarray.DataArray] or len(data.shape) < 2:
        raise ValueError(f'The data parameter in a  call to multislice_layout is'
                         f' not 2D or higher or is not an '
                         f'xarray.core.dataarray.DataArray object. data = '
                         f'{data}')

    if type(main_dims) in [list, tuple]:
        if len(main_dims) not in [2, 3]:  # If main axes is given check length
            raise ValueError(f'The main_dims parameter in a call to '
                             f'multislice_layout has the wrong number of '
                             f'elements. It expects a 2 or 3 element list/tuple'
                             f' but instead got {main_dims}')
        elif any([dim not in data.dims for dim in
                  main_dims]):  # elements of 'main_dims' are dimensions of data
            raise ValueError(f'In a call to multislice_layout one or more of '
                             f'main_dims is not a dimension of data. '
                             f'main_dims : {main_dims}, data dimensions : '
                             f'{data.dims}')
    else:
        main_dims = [dim for i, dim in enumerate(data.dims) if i < 3]

        # close figures to prevent to many being open if asked.
    if close_figures:
        plt.close()

    # generate the layout and return the required axes.
    fig, axes2D, axes1D, axesC = _generate_layout(
        main_dims=main_dims,
        extra_dims=[dim for dim in data.dims if dim not in main_dims],
        fig=fig)

    # Make the 2D plots
    plots = {}
    plot_parameters = {}
    cursors = {}
    for i, (axis_dims, axis) in enumerate(axes2D.items()):
        if i == 0 and len(main_dims) == 3:
            extent = [float(data.coords[main_dims[2]].max()),
                      float(data.coords[main_dims[2]].min()),
                      float(data.coords[main_dims[1]].min()),
                      float(data.coords[main_dims[1]].max())]
            transpose = False
            flip_dims = [main_dims[2]]
        elif i == 1 and len(main_dims) == 3:
            extent = [float(data.coords[main_dims[0]].min()),
                      float(data.coords[main_dims[0]].max()),
                      float(data.coords[main_dims[2]].max()),
                      float(data.coords[main_dims[2]].min())]
            transpose = True
            flip_dims = []
        elif i == 2 or len(main_dims) == 2:
            extent = [float(data.coords[main_dims[0]].min()),
                      float(data.coords[main_dims[0]].max()),
                      float(data.coords[main_dims[1]].min()),
                      float(data.coords[main_dims[1]].max())]
            transpose = True
            flip_dims = []

        _slice = data.sel({d: data.coords[d].mean()
                           for d in data.dims
                           if d not in axis_dims}, method='nearest')
        if transpose:
            _slice = _slice.T
        for flip_dim in flip_dims:
            _slice = _slice.isel({flip_dim: slice(None, None, -1)})

        plots[axis_dims] = axis.imshow(_slice, extent=extent, aspect='auto')
        plot_parameters[axis_dims] = {'transpose': transpose,
                                      'flip_dims': flip_dims}
        cursors[axis_dims] = {
            'x': [axis.axvline(
                x=data.coords[axis_dims[0]].mean(),
                color=cursor_colours[main_dims.index(axis_dims[1])],
                linestyle='--')],
            'y': [axis.axhline(
                y=data.coords[axis_dims[1]].mean(),
                color=cursor_colours[main_dims.index(axis_dims[0])],
                linestyle='--')]}

    # Make the 1Dplots
    for i, (axis_dims, axis) in enumerate(axes1D.items()):
        plots[axis_dims] = axis.plot(
            data.coords[axis_dims[0]],
            data.sel({d: data.coords[d].mean()
                      for d in data.dims
                      if d not in axis_dims}, method="nearest"),
            color=cursor_colours[i])[0]
        plot_parameters[axis_dims] = {'transpose': False, 'flip_dims': []}
        cursors[axis_dims] = {
            'x': [axis.axvline(x=data.coords[axis_dims[0]].mean(),
                               color='gray', linestyle='--')],
            'y': []}
    # Make the controls
    for i, dim in enumerate(data.dims):
        # create the slider
        slider = LinkableSlider(axesC[dim]['slider'], plots=plots,
                                cursors=cursors,
                                plot_parameters=plot_parameters,
                                label=dim, valmin=float(data.coords[dim].min()),
                                valmax=float(data.coords[dim].max()),
                                valstep=data.coords[dim],
                                valinit=float(data.coords[dim].mean()),
                                facecolor='lightgrey')
        slider.label.set_size(fontsize * fig_scale)
        slider.valtext.set_size(fontsize * fig_scale)
        slider.on_changed(_slider_update_factory(slider, dim))

        # Create the buttons
        left_button = mpl.widgets.Button(axesC[dim]['left_button'], '<',
                                         color='lightgray')
        left_button.on_clicked(_button_update_factory(slider,
                                                      direction='backward'))
        right_button = mpl.widgets.Button(axesC[dim]['right_button'], '>',
                                          color='lightgray')
        right_button.on_clicked(_button_update_factory(slider,
                                                       direction='forward'))

        if dim in controls.keys():  # create links to any existing sliders
            for link_slider in controls[dim]:
                link_slider.linked_sliders.append(slider)
                slider.linked_sliders.append(link_slider)
            controls[dim].append(slider)
            buttons[dim]['left'].append(left_button)
            buttons[dim]['right'].append(right_button)
        else:
            controls[dim] = [slider]
            buttons[dim] = {'left': [left_button], 'right': [right_button]}
    # Need to return 'controls' and 'buttons' otherwise the sliders and buttons
    # are inactive.
    return fig, controls, buttons
