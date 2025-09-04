"""
general submodule

This module contains general functions for customizing Matplotlib plots.

Functions:
    - check_font_size: Check the font size of labels and ticks on a plot.
    - check_axes_size: Check the size of axes in a plot.
    - check_all_rcParams: Check all Matplotlib rcParams settings.
    - get_ploted_figure: Get the current figure and axis objects.
    - get_points_on_markers_boundary: Get points on the boundary of markers in a plot.
    - add_arrow: Add an arrow with optional text annotation to a plot.
    - add_color_band: Add semi-transparent color band to axes.
    - add_circle_number: Add a numbered circle to a plot.
    - general_modify_band_plot: Modify the appearance of a band structure plot.
    - generate_gradient_colors: Generate a list of gradient colors.
    - general_modify_line: Modify the appearance of lines in a plot.
    - general_add_vlines_hlines: Add vertical and horizontal lines to a plot.
    - general_modify_legend: Modify the appearance of a legend.
    - general_margin_bin: Adjust the margins and binning for the x and y axes.
    - general_font: Customize global plotting settings for Matplotlib.
    - general_axes: Modifies axis properties like ticks, labels, and limits.
    - general_subplots_adjust: Adjusts the figure's subplots, size, and spacing.
    - general_adjust_text: Adjusts text labels to avoid overlap in a plot.
    - general_set_all_rcParams: Set multiple Matplotlib rcParams at once.
"""


from ase import Atoms

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import to_rgb
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, MaxNLocator
from matplotlib.patches import Ellipse

from brokenaxes import brokenaxes

from typing import List, Tuple, Union

from itertools import cycle
    
from adjustText import adjust_text

def check_font_size(ax: plt.Axes):
    x_label_fontsize = ax.xaxis.get_label().get_fontsize()
    y_label_fontsize = ax.yaxis.get_label().get_fontsize()
    print(f"Font size of x-axis label: {x_label_fontsize}")
    print(f"Font size of y-axis label: {y_label_fontsize}")

    x_tick_labels = ax.get_xticklabels()
    y_tick_labels = ax.get_yticklabels()

    if x_tick_labels:
        x_tick_font = x_tick_labels[0].get_fontproperties()
        print(f"Font name of the first major tick label on the x-axis: {x_tick_font.get_name()}")
        print(f"Font size of the first major tick label on the x-axis: {x_tick_font.get_size()}")

    if y_tick_labels:
        y_tick_font = y_tick_labels[0].get_fontproperties()
        print(f"Font name of the first major tick label on the y-axis: {y_tick_font.get_name()}")
        print(f"Font size of the first major tick label on the y-axis: {y_tick_font.get_size()}")

def check_axes_size(fig, ax1):
     # 1. 获取 ax1 在 Figure 中的位置和大小（归一化坐标，单位 0~1）
    bbox = ax1.get_position()
    print("ax1 in figure (relative):", bbox)

    # 2. 获取 figure 的物理尺寸（英寸）
    fig_w, fig_h = fig.get_size_inches()

    # 3. 计算 ax1 的实际尺寸（英寸）
    ax1_w_inch = bbox.width * fig_w
    ax1_h_inch = bbox.height * fig_h
    print(f"ax1 size in inches: {ax1_w_inch:.3f} x {ax1_h_inch:.3f}")

    # 4. 如果需要换算为像素（可选）
    dpi = fig.dpi
    print(f"ax1 size in pixels: {ax1_w_inch * dpi:.1f} x {ax1_h_inch * dpi:.1f}")

def check_all_rcParams():
    for key, value in plt.rcParams.items():
        print(f"{key:40} = {value}")

def get_ploted_figure():
    fig = plt.gcf()
    ax = plt.gca()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    return fig, ax, xlim, ylim

def get_points_on_markers_boundary(ax: plt.Axes = None, x: list = None, y: list = None, markersize_pts: float =20, ratio: float =1.0, num_points: int=20) -> np.ndarray:
    """Generates points around marker boundaries to avoid text overlap.

    Calculates a circular boundary around each marker point to help prevent
    text labels from overlapping with plot markers.

    Args:
        ax: Matplotlib axes object.
        x: X-coordinates of markers.
        y: Y-coordinates of markers.
        markersize_pts: Marker size in points.
        ratio: Scaling factor for marker boundary.
        num_points: Number of points to generate around each marker.

    Returns:
        np.ndarray: Array of (x,y) points around all marker boundaries.
    """
    fig = ax.figure
    ax_pos = ax.get_position()
    fig_width, fig_height = fig.get_size_inches()
    ax_width_inch = ax_pos.width * fig_width
    ax_height_inch = ax_pos.height * fig_height

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xlength = x1 - x0
    ylength = y1 - y0
    
    # radius, so /2
    markersize_inch = markersize_pts / 72 / 2 
    markersize_xlength = xlength * markersize_inch / ax_width_inch
    markersize_ylength = ylength * markersize_inch / ax_height_inch 
    
    angles = np.linspace(0, 2*np.pi, num_points)
    avoid_points = []
    for xi, yi in zip(x, y):
        px = xi + markersize_xlength * np.cos(angles)
        py = yi + markersize_ylength * np.sin(angles)
        avoid_points.extend(zip(px, py))
    
    return np.array(avoid_points)

def add_color_band(ax: Axes = None, extent: list = None, gradient: list = None, cmap: str = 'coolwarm', alpha: float = 0.5, origin: str='lower') -> None:
    """Adds a colored gradient band to an existing matplotlib Axes object.

    This function creates a semi-transparent color band on the specified axes using
    the given gradient colors or colormap. The band's position and size are determined
    by the extent parameter.

    Args:
        ax (Axes, optional): Matplotlib Axes object to add the color band to. 
            If None, uses current axes. Defaults to None.
        extent (list, optional): List of [xmin, xmax, ymin, ymax] specifying the 
            bounds of the color band in data coordinates. Defaults to None.
        gradient (list, optional): Array of color values for the gradient. 
            If None, uses the specified cmap. Defaults to None.
        cmap (str, optional): Matplotlib colormap name to use if gradient is None. 
            Defaults to 'coolwarm'.
        alpha (float, optional): Transparency of the color band (0-1). 
            Defaults to 0.5.
        origin (str, optional): Origin of the gradient image. Matches matplotlib's 
            imshow origin parameter. Defaults to 'lower'.

    Returns:
        None: The function modifies the Axes object in place.

    Example:
        >>> import matplotlib.pyplot as plt
        >>> fig, ax = plt.subplots()
        >>> ax.plot([0, 1], [0, 1])  # Some example plot
        >>> xmin, xmax = ax.get_xlim()
        >>> ymin, ymax = ax.get_ylim()
        >>> extent = [1/13, 1/11, ymin, ymax]
        >>> gradient_colors = generate_gradient_colors(if_cmap_color=True, 
        ...                                          if_reshape=True, 
        ...                                          reshape_M_N_L=[1, 1000, 4], 
        ...                                          if_reverse=True)
        >>> add_color_band(ax=ax, extent=extent, gradient=gradient_colors, alpha=0.5)
        >>> plt.show()

    Note:
        The function requires either a predefined gradient array or a valid cmap name.
        The gradient array should be shaped appropriately for matplotlib's imshow.
    """

    # Fill the gradient region
    ax.imshow(X=gradient, extent=extent, aspect='auto',
            cmap=cmap, origin=origin, alpha=alpha)

def add_circle_number(ax: Axes = None, positions: List[float]=None, color: str='steelblue', number: int=0,
                               radiusx_ratio: float = 0.05, 
                               markersize_pts: float = None, 
                               text_y_offset_inch: float = -0.0267, lw: float = 2, fontsize: int = 24,
                               **circle_kwargs) -> None:
    
    """Draws a numbered hollow ellipse at a specific data coordinate on a matplotlib axis.

    This function compensates for figure and axes aspect ratio distortions and allows
    precise placement of text annotations (numbers) inside ellipse markers.

    This function must be used after fixing the figure and axes size

    Args:
        positions (List[float]): A list of two elements [x, y] specifying the center of the ellipse.
        ax (Axes): A matplotlib Axes object to draw on.
        color (str): Color of the ellipse edge and number text. Defaults to 'steelblue'.
        number (int): The number to display inside the ellipse. Defaults to 0.
        radiusx_ratio (float): Horizontal radius of ellipse as a fraction of x-axis length. Defaults to 0.05.
        markersize_pts (float): Size of the marker in points. If None, uses radiusx_ratio to calculate size. Defaults to None.
        text_y_offset_inch (float): Vertical offset of the number text in inches. Defaults to -0.0267.
        lw (float): Line width of the ellipse edge. Defaults to 2.
        fontsize (int): Font size of the number text. Defaults to 24.
        **circle_kwargs: Additional keyword arguments passed to the Ellipse patch.

    Raises:
        ValueError: If `positions` is None, `ax` is None, or `positions` does not contain exactly two values.

    Examples:
        Basic usage:
        >>> fig, ax = plt.subplots()
        >>> add_circle_number([0.5, 0.5], ax, number=1, color='red')

        Annotating energy diagram features:
        >>> add_circle_number(
        ...     [max_fcc_layer2-0.01, min_energy_diff2+0.5*(max_energy_diff2-min_energy_diff2)],
        ...     ax, color=color4, radiusx_ratio=0.017, lw=2, number=1)
        >>> add_circle_number(
        ...     [max_fcc_layer2-0.01, min_energy_diff+0.5*(min_energy_diff2-min_energy_diff)],
        ...     ax, color=color4, radiusx_ratio=0.017, lw=2, number=2)

    Note:
        The function automatically compensates for:
        - Axes aspect ratio distortions
        - Figure-to-axes size scaling
        - Data coordinate transformations
        - Font rendering offsets
    """
    if positions is None:
        raise ValueError("positions must be provided")
    if ax is None:
        raise ValueError("ax must be provided")
    if len(positions) != 2:
        raise ValueError("positions must be a list of two elements [x, y]")

    ax_pos = ax.get_position()
    fig, _, (x0,x1), (y0,y1) = get_ploted_figure()
    
    fig_width, fig_height = fig.get_size_inches()
    ax_width_inches = ax_pos.width * fig_width
    ax_height_inches = ax_pos.height * fig_height
    xlength = x1 - x0
    ylength = y1 - y0

    aspect_ratio = ax_height_inches / ax_width_inches
    text_y_offset = text_y_offset_inch / ax_height_inches * ylength

    # For plot, tranform to data coordinates
    if markersize_pts is None:
        radiusx = radiusx_ratio * xlength
        # First trnasform to inch unit, then to data coordinates
        radiusy = radiusx_ratio / aspect_ratio * ylength
    else:
        # radius, so /2
        markersize_inch = markersize_pts / 72  / 2 
        radiusx = xlength * markersize_inch / ax_width_inches
        radiusy = ylength * markersize_inch / ax_height_inches 
        
    ax.add_patch(Ellipse(positions, width=2*radiusx, height=2*radiusy, edgecolor=color, facecolor='none', lw=lw,
                         **circle_kwargs))
    ax.text(positions[0], positions[1]+text_y_offset, str(number), fontsize=fontsize, ha='center', va='center', color=color)

def add_arrow(ax: plt.Axes = None, text: str = '', start: list = None, end: list = None,
              xycoords: str = 'data', textcoords: str = None, 
              arrowprops: dict = dict(arrowstyle="simple", color='red'), 
              zorder: int = 0, **kwargs) -> plt.Annotation:
    """
    Add an arrow between two points on a matplotlib Axes, optionally with text.

    This function is a wrapper around Matplotlib's annotation mechanism,
    specialized for drawing arrows from a starting point to an ending point.
    A text label will be added at the starting point.

    Args:
        ax (plt.Axes): The target Matplotlib Axes.
        text (str, optional): Annotation text to display near the arrow. Defaults to '', draw a single arrow.
        start (list[float], optional): Starting point of the arrow (x, y). Defaults to None. 
        end (list[float], optional): Ending point of the arrow (x, y). Defaults to None.
        xycoords (str or tuple or callable, optional): Coordinate system for 
            `start` and `end`. Defaults to `'data'`. Supported values include:

            - 'figure points': Points from the lower left of the figure.
            - 'figure pixels': Pixels from the lower left of the figure.
            - 'figure fraction': Fraction of the figure from the lower left.
            - 'subfigure points': Points from the lower left of the subfigure.
            - 'subfigure pixels': Pixels from the lower left of the subfigure.
            - 'subfigure fraction': Fraction of the subfigure from the lower left.
            - 'axes points': Points from the lower left corner of the Axes.
            - 'axes pixels': Pixels from the lower left corner of the Axes.
            - 'axes fraction': Fraction of the Axes from the lower left.
            - 'data': Data coordinates of the object being annotated 
              (default).
            - 'polar': (theta, r) if not in native 'data' coordinates.

            Note:
                'subfigure pixels' and 'figure pixels' are equivalent for the 
                parent figure. To ensure compatibility when using subfigures, 
                prefer 'subfigure pixels'.
        textcoords (str or tuple or callable, optional): Coordinate system for 
            the text position. If None, defaults to the value of `xycoords`.
            All `xycoords` values are valid, plus the following:

            - 'offset points': Offset in points from the xy value.
            - 'offset pixels': Offset in pixels from the xy value.
            - 'offset fontsize': Offset relative to the current font size 
              from the xy value.
        arrowprops (dict, optional): Properties used to draw a 
            .FancyArrowPatch arrow between the start and end positions. 
            Defaults to None, i.e., no arrow is drawn.

            There are two ways to specify arrows:

            Simple arrow (no 'arrowstyle' key):

            - width (float): Width of the arrow in points.
            - headwidth (float): Width of the base of the arrow head in points.
            - headlength (float): Length of the arrow head in points.
            - shrink (float): Fraction of the total arrow length to shrink from both ends.
            - ?: Any valid .FancyArrowPatch property.

            The arrow is attached to the edge of the text box. Its exact 
            position (corner or center) depends on where it points.

            Fancy arrow (when 'arrowstyle' is provided):

            - arrowstyle (str): Arrow style.
            - connectionstyle (str): Connection style.
            - relpos (tuple): Relative coordinates of the text box where 
              the arrow starts. Defaults to (0.5, 0.5).
            - patchA: Defaults to bounding box of the text.
            - patchB: Defaults to None.
            - shrinkA (float): Shrink factor in points. Defaults to 2.
            - shrinkB (float): Shrink factor in points. Defaults to 2.
            - mutation_scale (float): Scaling of mutation, defaults to text size (points).
            - mutation_aspect (float): Default is 1.
            - ?: Any valid .FancyArrowPatch property.

            Note:
                The starting point of the arrow is defined by relpos, where 
                (0, 0) is the lower left corner of the text box and (1, 1) is 
                the upper right. Values outside [0, 1] are supported and 
                indicate points outside the text box. Default is (0.5, 0.5), 
                i.e., centered.
        zorder (int, optional): Drawing order of the arrow/text. Defaults to 0.
        **kwargs: Additional keyword arguments passed to `.Text`.

    Returns:
        matplotlib.text.Annotation: The created annotation object.

    Example:
    ```python
        import matplotlib.pyplot as plt

        fig, ax = my_plot()

        # Add an arrow from (0.2, 0.2) to (0.8, 0.8) with high z-order
        add_arrow(ax=ax, start=[0.2, 0.2], end=[0.8, 0.8], zorder=10)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        plt.show()
    ```
    
    See Also:
        matplotlib.axes.Axes.annotate : Underlying annotation function.
        matplotlib.patches.FancyArrowPatch : Arrow drawing object.
    """
    temp = ax.annotate(text, xy=end, xytext=start, 
                xycoords=xycoords, textcoords=textcoords,
                arrowprops=arrowprops, zorder=zorder, **kwargs)
    return temp

def general_adjust_text(texts: list = None, ax=None, x: list = None, y: list = None, 
                        ensure_inside_axes: bool = True, iter_lim: int = 1000, time_lim: float = None,
                        expand: tuple = (1.0, 1.0), force_static: tuple = (0.1, 0.2), force_text: tuple = (0.1, 0.2),
                        if_strict_ensure_inside_axes: bool = True, density: int = 50,
                        if_avoid_markers: bool = True, markersize: float = 20.0, ratio_markers = 1.0, num_points_per_marker: int = 20,
                        **kwargs
                        ):
    """Adjusts text positions to avoid overlaps with markers and axes boundaries.

    Uses adjust_text algorithm to automatically reposition text labels to avoid collisions with
    specified points, plot markers, and axes boundaries.

    Args:
        texts: List of matplotlib.text.Text objects to adjust.
        ax: Matplotlib axes object.
        x: X-coordinates of points to avoid.
        y: Y-coordinates of points to avoid.
        ensure_inside_axes: Keep text within axes boundaries.
        iter_lim: Maximum iterations for adjustment algorithm.
        time_lim: Maximum time (seconds) for adjustment.
        expand: Expansion factors for text bounding boxes.
        force_static: Repulsion force from static points.
        force_text: Repulsion force between text labels.
        if_strict_ensure_inside_axes: Add boundary points to avoid.
        density: Number of boundary points to generate.
        if_avoid_markers: Avoid overlap with plot markers.
        markersize: Size of markers to avoid.
        ratio_markers: Ratio of marker size to avoid zone.
        num_points_per_marker: Points to generate around each marker.
        **kwargs: Additional arguments passed to adjust_text.

    Raises:
        ValueError: If x and y lengths don't match, or if either is None.

    Note:
        Requires the adjustText package (https://github.com/Phlya/adjustText).

    Example:
        ```python
        from mymetal.universal.plot.general import general_adjust_text
        import numpy as np
        from mymetal.universal.plot.plot import my_plot

        # Typical usage in a plotting function to avoid label overlaps
        fig, ax = my_plot() 
        texts = []
        positions = np.random.rand(10, 2)  # Example data points
        ax.plot(positions[:, 0], positions[:, 1], 'o')  

        # Create some text labels
        for i, (x, y) in enumerate(positions):
            texts.append(ax.text(x, y, f'Label {i+1}', ha='center', va='center'))
        
        # Adjust text positions to avoid overlaps
        general_adjust_text(
            texts=texts,
            ax=ax,
            x=positions[:, 0],
            y=positions[:, 1],
            ensure_inside_axes=True,
            if_avoid_markers=True,
            # Optional parameters with defaults:
            # expand=(1.0, 1.0),
            # force_static=(0.1, 0.2),
            # force_text=(0.1, 0.2),
            # if_strict_ensure_inside_axes=True,
            # density=50,
            # markersize=20.0,
            # ratio_markers=1.0,
            # num_points_per_marker=20
        )
        
        ```

    """
    # markers & check
    if x is not None and y is not None:
        if len(x) != len(y):
            raise ValueError('x and y must have the same length.')
        else:
            avoid_points = np.column_stack([x, y])
    if x is None or y is None:
        raise ValueError('Both x and y must be provided for avoid points.')

    # boundary points
    # TODO: ratio to boundary
    avoid_boundary_points = []
    if if_strict_ensure_inside_axes:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        top = np.column_stack([np.linspace(xmin, xmax, density), 
                            np.full(density, ymax)])
        bottom = np.column_stack([np.linspace(xmin, xmax, density), 
                                np.full(density, ymin)])
        left = np.column_stack([np.full(density, xmin),
                            np.linspace(ymin, ymax, density)])
        right = np.column_stack([np.full(density, xmax),
                            np.linspace(ymin, ymax, density)])
        avoid_boundary_points = np.vstack([top, bottom, left, right])

    avoid_marker_points = []
    if if_avoid_markers:
        avoid_marker_points = get_points_on_markers_boundary(ax, x, y, markersize_pts=markersize, ratio=ratio_markers, num_points=num_points_per_marker)

    avoid_full_points = np.vstack([avoid_points, avoid_boundary_points, avoid_marker_points])
    avoid_full_x = avoid_full_points[:, 0]
    avoid_full_y = avoid_full_points[:, 1]
    adjust_text(
        texts = texts,
        ax=ax,
        x = avoid_full_x,
        y = avoid_full_y,
        expand = expand,
        force_static=force_static,
        force_text=force_text,
        ensure_inside_axes=ensure_inside_axes,
        iter_lim= iter_lim,  
        time_lim=time_lim,
        **kwargs)

def general_modify_band_plot(ax: plt.Axes) -> plt.Axes:
    xticks = ax.get_xticks()
    xticklabels = ax.get_xticklabels()
    xlabels = ax.get_xlabel()
    ax.tick_params(axis='x', which='both', bottom=False)
    general_add_vlines_hlines(ax, xticks)
    return ax

def generate_gradient_colors(start_color: str = 'red', end_color: str = 'blue', num_colors: int = 2,
                             if_cmap_color: bool = False, cmap_color: str = 'coolwarm', if_reshape: bool = False, reshape_M_N_L: list = [1, 1000, 4],
                             if_reverse: bool = False) -> list:
    """Generates a list of colors forming a gradient between two colors or using a colormap.

    This function can create either a simple linear gradient between two colors or
    sample colors from a matplotlib colormap. The output can be reshaped into a
    multi-dimensional array and optionally reversed.

    Args:
        start_color (str, optional): Starting color for the gradient (named color or hex).
            Defaults to 'red'.
        end_color (str, optional): Ending color for the gradient (named color or hex).
            Defaults to 'blue'.
        num_colors (int, optional): Number of colors to generate in the gradient.
            Defaults to 2.
        if_cmap_color (bool, optional): If True, uses a colormap instead of start/end colors.
            Defaults to False.
        cmap_color (str, optional): Name of matplotlib colormap to use if if_cmap_color is True.
            Defaults to 'coolwarm'.
        if_reshape (bool, optional): If True, reshapes output array according to reshape_M_N_L.
            Defaults to False.
        reshape_M_N_L (list, optional): Dimensions for reshaping output as [M, N, L].
            Defaults to [1, 1000, 4].
        if_reverse (bool, optional): If True, reverses the order of colors in the gradient.
            Defaults to False.

    Returns:
        list: Array of RGBA color values. Shape depends on if_reshape parameter:
            - If not reshaped: (num_colors, 4)
            - If reshaped: (M, N, L) as specified

    Example:
        >>> # Simple two-color gradient
        >>> colors = generate_gradient_colors('red', 'blue', 5)
        >>> 
        >>> # Using a colormap with 100 colors
        >>> colors = generate_gradient_colors(if_cmap_color=True, cmap_color='viridis', num_colors=100)
        >>>
        >>> # Reshaped gradient for use with imshow
        >>> gradient = generate_gradient_colors(
        ...     if_cmap_color=True,
        ...     cmap_color='coolwarm',
        ...     if_reshape=True,
        ...     reshape_M_N_L=[1, 100, 4],
        ...     if_reverse=True
        ... )
        >>> plt.imshow(gradient, aspect='auto')
        >>> plt.show()

    Note:
        When if_reshape=True, the product of M*N must equal num_colors (or the default
        num_colors value if if_reshape is True but num_colors wasn't specified).
        The L dimension should typically be 4 (for RGBA values).
    """
    if if_reshape:
        M = reshape_M_N_L[0]
        N = reshape_M_N_L[1]
        L = reshape_M_N_L[2]
        num_colors = M * N
    if if_cmap_color:
        cmap = plt.get_cmap(cmap_color)
        values = np.linspace(0, 1, num_colors)
        color_list = cmap(values)
    else:
        cmap = mcolors.LinearSegmentedColormap.from_list("gradient", [start_color, end_color])
        color_list = [cmap(i / (num_colors - 1)) for i in range(num_colors)]
    if if_reverse:
        temp_color_list = color_list.copy()
        color_list = temp_color_list[::-1]
    if if_reshape:
        temp_color_list = np.array(color_list.copy())
        color_list = temp_color_list.reshape(M, N, L)
    return color_list

def general_modify_line(ax: plt.Axes,
                    remove_dashed_line: bool = True,
                    remove_line_style: list = ['--', ':'],
                    color_list: list = ['black', 'red', 'blue', 'green', 'orange', 'purple', 
                                        'brown', 'pink', 'gray', 'olive', 'cyan', 'lime', 
                                        'teal', 'lavender', 'tan', 'salmon', 'gold', 'lightcoral', 
                                        'skyblue', 'darkblue', 'darkred', 'darkgreen', 'darkorange', 
                                        ],
                    if_gradient_color: bool = False,
                    gradient_colors: list=['red', 'blue'],
                    if_cmap_color: bool = False,
                    cmap_color: str = 'coolwarm',
                    if_change_color: bool = False,
                    ):
    all_lines = ax.get_lines()
    all_index = [index for index, line in enumerate(all_lines)]
    removed_index = []
    after_removed_lines = []
    if remove_dashed_line:
        for index, line in enumerate(all_lines):
            if line.get_linestyle() in remove_line_style:
                line.set_visible(False)
                removed_index.append(index)
            else:
                after_removed_lines.append(line)
    if if_gradient_color:
        color_list = generate_gradient_colors(gradient_colors[0], gradient_colors[1], len(after_removed_lines),
                                              if_cmap_color, cmap_color)
    color_cycle = cycle(color_list)
    if if_change_color:
        for index, line in enumerate(after_removed_lines):
            line.set_color(next(color_cycle))
    return ax

def general_add_vlines_hlines(  ax: plt.Axes,
                                vlines: list = [],
                                hlines: list = [],
                                vlines_colors: list = ['gray'],
                                hlines_colors: list = ['gray'],
                                zorder_vlines: list = [-1],
                                zorder_hlines: list = [-1],
                                linestyle: str='--'):
    # add vlines, hlines
    color_cycle = cycle(vlines_colors)
    zorder_cycles=cycle(zorder_vlines)
    for vline in vlines:
        ax.axvline(x=vline, color=next(color_cycle), zorder=next(zorder_cycles), linestyle= linestyle)
    color_cycle = cycle(hlines_colors)
    zorder_cycles=cycle(zorder_hlines)
    for hline in hlines:
        ax.axhline(y=hline, color=next(color_cycle), zorder=next(zorder_cycles),  linestyle= linestyle)

def general_modify_legend(legend, boxstyle: str = 'Square',
                        pad: float = 0.5, linewidth: float = 2.5,
                        edgecolor: str = 'black', facecolor: str = 'white', alpha: float = 1):
    """
    Customizes the appearance of a legend.

    Args:
        legend: The legend object to modify.

    Returns:
        None

    Example:
            >>> from mymetal.universal.plot.plot import my_plot
            >>> from mymetal.universal.plot.general import general_modify_legend
            >>> fig, axes = my_plot(left = 1.6, grid=False, fig_sharex=False)
            >>> handles_left, labels_left = ax.get_legend_handles_labels()
            >>> handles_right, labels_right = ax2.get_legend_handles_labels()
            >>> handles = handles_left + handles_right
            >>> labels = labels_left + labels_right
            >>> legend = ax.legend(handles, labels, loc='upper right',
            ...                    bbox_to_anchor=(0.95, 0.7),
            ...                    columnspacing=0.0, ncol=1)
            >>> general_modify_legend(ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.05), fontsize = 24))
    """
    legend.get_frame().set_boxstyle(boxstyle, pad=pad)
    legend.get_frame().set_linewidth(linewidth)
    legend.get_frame().set_edgecolor(edgecolor)
    legend.get_frame().set_facecolor(facecolor)
    legend.get_frame().set_alpha(alpha)
    
def general_margin_bin(axis,
                        y_margin = 0.1, 
                        y_nbins = 3,
                        x_margin = 0.1, 
                        x_nbins = 4, 
                        prune = 'both',
                        if_autoscale = True):
    """
    Adjusts the margins and binning for the x and y axes.

    Args:
        axis: The axis object to modify.
        y_margin (float): The y-axis margin. Default is 0.1.
        y_nbins (int): Number of bins for the y-axis. Default is 3.
        x_margin (float): The x-axis margin. Default is 0.1.
        x_nbins (int): Number of bins for the x-axis. Default is 4.
        prune (str): Prune setting for axis ticks. Default is 'both'.
        if_autoscale (bool): Whether to autoscale the axis. Default is True.

    Returns:
        None
    """
    axis.margins(x=x_margin, y=y_margin)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins, prune=prune))
    axis.xaxis.set_major_locator(MaxNLocator(nbins=x_nbins, prune=prune))
    if if_autoscale:
        # Force axis limits to update based on the new margins
        axis.autoscale(enable=True, axis='both', tight=False)

def general_font( grid: bool = True, grid_linewidth: float = 2, grid_linestyle = '--',
                 if_ploted: bool = False,
                  fontsize: int = 28, markersize: int = 20, linewidth: int = 3, 
                  markeredgewidth: int = 3, legend_fontsize: int = 24, 
                  markerfacecolor: str = 'white', fontname: str = 'Arial',
                  loc: str = 'upper right', frameon: bool = True, 
                  borderpad: float=0.0, labelspacing: float=0.5, columnspacing:float=0.5,
                  constrained_layout: bool = False) -> None:
    """
    Customize global plotting settings for Matplotlib, including font styles, line properties, 
    grid appearance, and legend settings. Optionally, apply these settings to an existing plot.

    Args:
        grid (bool, optional): Whether to display grid lines on the plot. Default is True.
        grid_linewidth (float, optional): Line width for grid lines. Default is 0.5.
        grid_linestyle (str, optional): Style of grid lines (e.g., '--', '-'). Default is '--'.
        if_ploted (bool, optional): Whether to apply the settings to the current plot if one exists. Default is False.
        fontsize (int, optional): Font size for text elements such as labels, title, and ticks. Default is 28.
        markersize (int, optional): Size of markers in the plot. Default is 20.
        linewidth (int, optional): Width of plot lines. Default is 3.
        markeredgewidth (int, optional): Width of marker edges. Default is 3.
        legend_fontsize (int, optional): Font size for the legend text. Default is 24.
        markerfacecolor (str, optional): Color of the marker face. Default is 'white'.
        fontname (str, optional): Name of the font to use for text elements. Default is 'Arial'.
        loc (str, optional): Location of the legend on the plot. Default is 'upper right'.
        frameon (bool, optional): Whether to display a frame around the legend. Default is True.
        borderpad (float, optional): Padding around the legend. Default is 0.0.
        labelspacing (float, optional): Spacing between legend labels. Default is 0.5.
        columnspacing (float, optional): Spacing between legend columns. Default is 0.5.

    Returns:
        None: This function modifies the global plotting settings and optionally the current plot.

    Example:
        >>> general_font(fontsize=20, legend_fontsize=18, if_ploted=True)
    """
    plt.rcParams['font.family'] = fontname                 #
    plt.rcParams['font.size'] = fontsize                   #
    plt.rcParams['axes.linewidth'] = linewidth             #
    plt.rcParams['axes.grid'] = grid                       #
    plt.rcParams['grid.linestyle'] = grid_linestyle        #
    plt.rcParams['grid.linewidth'] = grid_linewidth        #
    plt.rcParams['lines.linewidth'] = linewidth            #
    plt.rcParams['lines.markersize'] = markersize          #
    plt.rcParams['lines.markeredgewidth'] = markeredgewidth #
    plt.rcParams['lines.markerfacecolor'] = markerfacecolor #
    
    # layout
    plt.rcParams['figure.constrained_layout.use'] = constrained_layout

    # legend     will done, because the legend will be re-drawed
    plt.rcParams['legend.loc'] = loc
    plt.rcParams['legend.fontsize'] = legend_fontsize
    plt.rcParams['legend.frameon'] = frameon
    plt.rcParams['legend.borderpad'] = borderpad
    plt.rcParams['legend.labelspacing'] = labelspacing
    plt.rcParams['legend.columnspacing'] = columnspacing

    # Math
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'  # 设置数学字体为 Arial

    # Savefig
    plt.rcParams['savefig.transparent'] = False
    plt.rcParams['savefig.facecolor'] = 'white'
    plt.rcParams['savefig.dpi'] = 300 

    # background
    plt.rcParams['axes.facecolor'] = 'white' 
    plt.rcParams['figure.facecolor'] = 'white'  

    if if_ploted:
        fig = plt.gcf()
        axis = plt.gca()
        if isinstance(axis, np.ndarray):
            axes = axis.flatten()
        else:  
            axes = [axis]
        for ax in axes:
            # change the font in tick_params
            plt.xticks(fontsize=fontsize, fontname= fontname)
            plt.yticks(fontsize=fontsize, fontname= fontname)
            # change grid
            ax.grid(linewidth=grid_linewidth, linestyle = grid_linestyle)
            ax.grid(grid)
            # re-generate title/xlabel/ylabel about fontsize, fontname
            ax.set_title(ax.get_title(), size=fontsize, fontname=fontname)
            ax.set_xlabel(ax.get_xlabel(), size=fontsize, fontname=fontname)
            ax.set_ylabel(ax.get_ylabel(), size=fontsize, fontname=fontname)
            # change the linewidth in ax.
            for temp in ['top', 'bottom', 'right', 'left']:
                ax.spines[temp].set_linewidth(linewidth) 
            # change line setting
            for line in ax.get_lines():
                line.set_linewidth(linewidth)
                line.set_markersize(markersize)
                line.set_markeredgewidth(markeredgewidth)
                line.set_markerfacecolor(markerfacecolor)
        
def general_axes(ax,
                labelpad: int = 15, 
                tick_pad: int = 10, 
                xlabel: str = 'X-axis Label',
                ylabel: str = 'Y-axis Label',
                xlim: Union[List[float], None] = None,
                ylim: Union[List[float], None] = None,
                if_set_lim: bool = True,
                if_close_right_top_tick: bool = False,
                if_ax2: bool = False,
                color: str = 'black',):
    """
    Modifies axis properties like ticks, labels, and limits.

    Args:
        ax: The axis object to modify.
        labelpad (int): Padding for axis labels. Default is 15.
        tick_pad (int): Padding for axis ticks. Default is 10.
        xlabel (str): Label for the x-axis. Default is 'X-axis Label'.
        ylabel (str): Label for the y-axis. Default is 'Y-axis Label'.
        xlim (Union[List[float], None]): x-axis limits. Default is None.
        ylim (Union[List[float], None]): y-axis limits. Default is None.

    Returns:
        None
    """
    if isinstance(ax, np.ndarray):
        axes = ax.flatten()
    else:  
        axes = [ax]

    # change the top/bottom/top/right linewidth in general_font() function.
    for axis in axes:
        axis.minorticks_on()
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
        axis.tick_params(which='major', direction='in', length=8, width=3.0, pad = tick_pad)
        axis.tick_params(which='minor', direction='in', length=4, width=3.0, pad = tick_pad)

        axis.set_xlabel(xlabel, labelpad=labelpad)
        axis.set_ylabel(ylabel, labelpad=labelpad)

        if if_close_right_top_tick:
            axis.tick_params(axis='y', which='both', right=False)
            axis.tick_params(axis='x', which='both', top=False)

        if if_set_lim:
            axis.set_xlim(xlim)
            axis.set_ylim(ylim)
        
        if if_ax2:
            axis.set_ylabel(ylabel, color=color, labelpad = labelpad)
            axis.tick_params(axis='y', labelcolor=color, direction='in', color = color)
            axis.spines['right'].set_color(color) 
            axis.tick_params(axis='y', which='both', color = color)

def general_subplots_adjust(fig: Figure,
                            one_fig_wh: List[float] = [10.72, 8.205], 
                            fig_subp: List[int] = [1, 1],   
                            axes_height: float = 5.89,
                            axes_width: float = 7.31,
                            left: float = 1.918, 
                            top: float = 0.9517):
    """
    Adjusts the figure's subplots, size, and spacing.

    Args:
        fig: The figure object to adjust.
        one_fig_wh (List[float]): Width and height of one figure. Default is [10.72, 8.205].
        fig_subp (List[int]): Number of subplots. Default is [1, 1].
        axes_height (float): Height of the axes. Default is 5.89.
        axes_width (float): Width of the axes. Default is 7.31.
        left (float): Left margin. Default is 1.918.
        top (float): Top margin. Default is 0.9517.

    Returns:
        None
    """
    fig_wh = [fig_subp[1]*one_fig_wh[0], fig_subp[0]*one_fig_wh[1]]
    # inch 10.72* 8.205 is one figsize
    left = left
    right = one_fig_wh[0]-left-axes_width
    top = top
    bottom = one_fig_wh[1]-top-axes_height
    wspace = (left+right)/axes_width  # axes_width 7.31 inch
    hspace = (top+bottom)/axes_height  # axes_height 5.89 inch

    fig.subplots_adjust(left=left/fig_wh[0], right=(fig_wh[0]-right)/fig_wh[0], 
                        top=(fig_wh[1]-top)/fig_wh[1], bottom=(bottom)/fig_wh[1],
                        hspace= (hspace), wspace=(wspace))

def general_set_all_rcParams(
        # axes
        axes_grid: bool = True,
        axes_labelpad: float = 15,
        axes_labelsize: float = 28,
        axes_linewidth: float = 3,
        axes_titlepad: float = 15,
        axes_xmargin: float = 0.1,
        axes_ymargin: float = 0.1,
        axes_zmargin: float = 0.1,
        # backend
        backend: str = 'module://matplotlib_inline.backend_inline',
        # boxplot, Todo

        # contour
        contour_linewidth: float = 3,

        # errorbar
        errorbar_capsize: float = 3,

        # figure
        figure_dpi: float = 300.0,
        figure_edgecolor: str = "white",
        figure_facecolor: str = "white",
        # control figsize by subp and one_figsize
        figure_subp: tuple = (1, 1),                # (rows, cols)
        figure_one_figsize: tuple = (10.72, 8.205), # inches, (width, height) 

        figure_frameon: bool = True,

        # control subplot by left, right, top, bottom, hspace, wspace
        figure_left: float = 1.918, # inches
        figure_top: float = 0.9517, # inches
        figure_wspace: float = 3.41, # inches
        figure_hspace: float = 2.315, # inches

        # font
        font_family = [
            "Arial",                # 无衬线，常用
            "Times New Roman",      # 有衬线，经典论文字体
            # "Computer Modern",      # LaTeX 默认字体
            # "Latin Modern Roman",   # LaTeX 常用衬线字体
            # "DejaVu Sans",          # Matplotlib 默认备选
            # "DejaVu Serif"          # Matplotlib 默认备选
        ],
        fontsize: float = 28,

        # grid
        grid_alpha: float = 1.0,
        grid_color: str = 'gray',
        grid_linestyle: str = '--',
        grid_linewidth: float = 2,

        # legend
        legend_loc: str = 'upper right', 
        legend_frameon: bool = True, 
        legend_borderpad: float=0.0,
        legend_labelspacing: float=0.5,
        legend_columnspacing:float=0.5,
        legend_edgecolor: str = 'black',
        legend_facecolor: str = 'white',
        legend_framealpha: float = 1.0,
        legend_boxstyle: str = 'Square',
        legend_pad: float = 0.5,
        legend_linewidth: float = 2.5,

        # lines
        lines_linestyle: str = '-',
        lines_linewidth: float = 3,
        lines_markeredgecolor = 'auto',
        lines_markeredgewidth: float = 3,
        lines_markerfacecolor = 'white',
        lines_markersize: float = 20,

        # mathtext
        mathtext_fontset: str = 'custom',
        mathtext_rm: str = 'Arial',  # 设置数学字体为 Arial
        mathtext_it: str = 'Arial:italic', # 设置数学斜体为 Arial Italic
        mathtext_bf: str = 'Arial:bold', # 设置数学粗体为 Arial Bold
        mathtext_sf: str = 'Arial', # 设置数学无衬线为 Arial
        mathtext_tt: str = 'Arial', # 设置数学等宽字体为 Arial
        mathtext_cal: str = 'Arial', # 设置数学手写体为 Arial
        mathtext_default: str = 'it', # it, rm, cal, sf, tt

        # patch
        patch_linewidth: float = 3.0,

        # savefig
        savefig_dpi: float = 300,
        savefig_format: str = 'png',
        savefig_transparent: bool = False,

        # x-axis tick settings
        xtick_bottom: bool = True,
        xtick_top: bool = False,
        xtick_direction: str = 'in',
        xtick_labelbottom: bool = True,
        xtick_labeltop: bool = False,
        xtick_labelsize: float = 28,
        xtick_major_bottom: bool = True,
        xtick_major_top: bool = True,
        xtick_major_size: float = 8,
        xtick_major_width: float = 3.0,
        xtick_major_pad: float = 10,
        xtick_minor_bottom: bool = True,
        xtick_minor_top: bool = True,
        xtick_minor_size: float = 4,
        xtick_minor_width: float = 3.0,
        xtick_minor_pad: float = 10,
        xtick_minor_ndivs: float = 2,
        xtick_minor_visible: bool = True,

        # y-axis tick settings
        ytick_left: bool = True,
        ytick_right: bool = False,
        ytick_direction: str = 'in',
        ytick_labelleft: bool = True,
        ytick_labelright: bool = False,
        ytick_labelsize: float = 28,
        ytick_major_left: bool = True,
        ytick_major_right: bool = True,
        ytick_major_size: float = 8,
        ytick_major_width: float = 3.0,
        ytick_major_pad: float = 10,
        ytick_minor_left: bool = True,
        ytick_minor_right: bool = True,
        ytick_minor_size: float = 4,
        ytick_minor_width: float = 3.0,
        ytick_minor_pad: float = 10,
        ytick_minor_ndivs: str = 2,
        ytick_minor_visible: bool = True,
):
    """
    Configure global matplotlib rcParams for consistent figure styling.

    This function can be used as a replacement for `my_plot()`, and the performance of the two 
    functions is essentially identical. It modifies nearly all rcParams that are typically adjusted 
    during plotting, but many parameters remain untouched.

    This function centralizes the control of matplotlib's default rendering style 
    for consistent and publication-quality figures.

    Example:
        >>> import matplotlib.pyplot as plt
        >>> from mymetal.universal.plot.general import general_all_rcParams
        >>> subp = (2, 2)
        >>> lg = general_all_rcParams(figure_subp=subp)
        >>> fig, axes = plt.subplots(subp[0], subp[1])
        >>> ax = axes.flatten()[0]
        >>> ax.plot([0, 1], [0, 1], marker='o', label='Line 1')
        >>> ax.plot([0, 1], [1, 0], marker='s', label='Line 2')
        >>> ax.set_xlabel("X Label")
        >>> ax.set_ylabel("Y Label")
        >>> rect = plt.Rectangle((0.1, 0.1), 0.5, 0.5)
        >>> ax.add_patch(rect)
        >>> lg(ax.legend())
        >>> ax.set_title('Demo of rcParams setup')
        >>> plt.savefig('demo_plot.png')

    Note:
        You need to call `plt._general_modify_legend(legend)` after creating a legend to apply
        the custom legend styling.
    """    

    # backend: inline, TkAgg, Agg, etc.
    plt.rcParams['axes.grid'] = axes_grid 
    plt.rcParams['axes.labelpad'] = axes_labelpad
    plt.rcParams['axes.labelsize'] = axes_labelsize
    plt.rcParams['axes.linewidth'] = axes_linewidth
    plt.rcParams['axes.titlepad'] = axes_titlepad

    plt.rcParams['axes.xmargin'] = axes_xmargin
    plt.rcParams['axes.ymargin'] = axes_ymargin
    plt.rcParams['axes.zmargin'] = axes_zmargin

    plt.rcParams['backend'] = backend

    plt.rcParams['contour.linewidth'] = contour_linewidth

    plt.rcParams['errorbar.capsize'] = errorbar_capsize
    
    # figure
    plt.rcParams['figure.dpi'] = figure_dpi
    plt.rcParams['figure.edgecolor'] = figure_edgecolor
    plt.rcParams['figure.facecolor'] = figure_facecolor

    # calculate figsize and subplot params
    figure_figsize = (figure_one_figsize[0] * figure_subp[1], figure_one_figsize[1] * figure_subp[0])
    plt.rcParams['figure.figsize'] = figure_figsize
    plt.rcParams['figure.frameon'] = figure_frameon
    figure_right = figure_wspace - figure_left
    figure_bottom = figure_hspace - figure_top
    figure_subplot_left = figure_left / figure_figsize[0]
    figure_subplot_right = 1 - figure_right / figure_figsize[0]
    figure_subplot_top  = 1 - figure_top / figure_figsize[1]
    figure_subplot_bottom = figure_bottom / figure_figsize[1]
    figure_one_axes_width = figure_one_figsize[0] - figure_wspace # inches
    figure_one_axes_height = figure_one_figsize[1] - figure_hspace # inches
    figure_subplot_wspace = figure_wspace / figure_one_axes_width 
    figure_subplot_hspace = figure_hspace / figure_one_axes_height
    plt.rcParams['figure.subplot.left']   = figure_subplot_left
    plt.rcParams['figure.subplot.right']  = figure_subplot_right
    plt.rcParams['figure.subplot.top']    = figure_subplot_top
    plt.rcParams['figure.subplot.bottom'] = figure_subplot_bottom
    plt.rcParams['figure.subplot.hspace'] = figure_subplot_hspace
    plt.rcParams['figure.subplot.wspace'] = figure_subplot_wspace

    # font
    plt.rcParams['font.family'] = font_family
    plt.rcParams['font.size'] = fontsize

    # grid
    plt.rcParams['grid.alpha'] = grid_alpha
    plt.rcParams['grid.color'] = grid_color
    plt.rcParams['grid.linestyle'] = grid_linestyle
    plt.rcParams['grid.linewidth'] = grid_linewidth

    # legend
    plt.rcParams['legend.loc'] = legend_loc
    plt.rcParams['legend.frameon'] = legend_frameon
    plt.rcParams['legend.borderpad'] = legend_borderpad
    plt.rcParams['legend.labelspacing'] = legend_labelspacing
    plt.rcParams['legend.columnspacing'] = legend_columnspacing
    plt.rcParams['legend.edgecolor'] = legend_edgecolor
    plt.rcParams['legend.facecolor'] = legend_facecolor
    plt.rcParams['legend.framealpha'] = legend_framealpha

    from mymetal.universal.plot.general import general_modify_legend
    def _general_modify_legend(legend, ):
        """general_modify_legend(ax.legend(loc='lower center', 
        bbox_to_anchor=(0.5, 0.05), fontsize = 24))
        """
        general_modify_legend(legend, boxstyle=legend_boxstyle, pad=legend_pad,
                            linewidth=legend_linewidth, 
                            edgecolor=legend_edgecolor,
                            facecolor=legend_facecolor,
                            alpha=legend_framealpha)
    
    
    # lines
    plt.rcParams['lines.linestyle'] = lines_linestyle
    plt.rcParams['lines.linewidth'] = lines_linewidth
    plt.rcParams['lines.markeredgecolor'] = lines_markeredgecolor
    plt.rcParams['lines.markeredgewidth'] = lines_markeredgewidth
    plt.rcParams['lines.markerfacecolor'] = lines_markerfacecolor
    plt.rcParams['lines.markersize'] = lines_markersize

    # mathtext
    plt.rcParams['mathtext.fontset'] = mathtext_fontset
    plt.rcParams['mathtext.rm'] = mathtext_rm
    plt.rcParams['mathtext.it'] = mathtext_it
    plt.rcParams['mathtext.bf'] = mathtext_bf
    plt.rcParams['mathtext.sf'] = mathtext_sf
    plt.rcParams['mathtext.tt'] = mathtext_tt
    plt.rcParams['mathtext.cal'] = mathtext_cal
    plt.rcParams['mathtext.default'] = mathtext_default

    # patch
    plt.rcParams['patch.linewidth'] = patch_linewidth

    # savefig
    plt.rcParams['savefig.dpi'] = savefig_dpi
    plt.rcParams['savefig.format'] = savefig_format
    plt.rcParams['savefig.transparent'] = savefig_transparent

    # x-axis tick settings
    plt.rcParams['xtick.bottom'] = xtick_bottom
    plt.rcParams['xtick.top'] = xtick_top
    plt.rcParams['xtick.direction'] = xtick_direction
    plt.rcParams['xtick.labelbottom'] = xtick_labelbottom
    plt.rcParams['xtick.labeltop'] = xtick_labeltop
    plt.rcParams['xtick.labelsize'] = xtick_labelsize
    plt.rcParams['xtick.major.bottom'] = xtick_major_bottom
    plt.rcParams['xtick.major.top'] = xtick_major_top
    plt.rcParams['xtick.major.size'] = xtick_major_size
    plt.rcParams['xtick.major.width'] = xtick_major_width
    plt.rcParams['xtick.major.pad'] = xtick_major_pad
    plt.rcParams['xtick.minor.bottom'] = xtick_minor_bottom
    plt.rcParams['xtick.minor.top'] = xtick_minor_top
    plt.rcParams['xtick.minor.size'] = xtick_minor_size
    plt.rcParams['xtick.minor.width'] = xtick_minor_width
    plt.rcParams['xtick.minor.pad'] = xtick_minor_pad
    plt.rcParams['xtick.minor.ndivs'] = xtick_minor_ndivs
    plt.rcParams['xtick.minor.visible'] = xtick_minor_visible

    # y-axis tick settings
    plt.rcParams['ytick.left'] = ytick_left
    plt.rcParams['ytick.right'] = ytick_right
    plt.rcParams['ytick.direction'] = ytick_direction
    plt.rcParams['ytick.labelleft'] = ytick_labelleft
    plt.rcParams['ytick.labelright'] = ytick_labelright
    plt.rcParams['ytick.labelsize'] = ytick_labelsize
    plt.rcParams['ytick.major.left'] = ytick_major_left
    plt.rcParams['ytick.major.right'] = ytick_major_right
    plt.rcParams['ytick.major.size'] = ytick_major_size
    plt.rcParams['ytick.major.width'] = ytick_major_width
    plt.rcParams['ytick.major.pad'] = ytick_major_pad
    plt.rcParams['ytick.minor.left'] = ytick_minor_left
    plt.rcParams['ytick.minor.right'] = ytick_minor_right
    plt.rcParams['ytick.minor.size'] = ytick_minor_size
    plt.rcParams['ytick.minor.width'] = ytick_minor_width
    plt.rcParams['ytick.minor.pad'] = ytick_minor_pad
    plt.rcParams['ytick.minor.ndivs'] = ytick_minor_ndivs
    plt.rcParams['ytick.minor.visible'] = ytick_minor_visible

    return _general_modify_legend