import matplotlib as mpl


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
