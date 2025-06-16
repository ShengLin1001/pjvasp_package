"""
render submodule

This module contains functions to render images from OVITO pipelines.

Functions:
    - my_render: Render an image from an OVITO pipeline.
"""


from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer, OSPRayRenderer
from ovito.modifiers import CommonNeighborAnalysisModifier, AssignColorModifier, ColorCodingModifier, ExpressionSelectionModifier, ClearSelectionModifier
from ovito.pipeline import Pipeline
import ovito


def my_render(pipeline=None, imagefile=None, size=(3200, 2400), renderer=TachyonRenderer(shadows=True, direct_light_intensity=1.1),
                      camera_dir = (1, 1, 1), viewtype = Viewport.Type.Ortho):
    """
    Renders an image from an OVITO pipeline using a specified camera direction and renderer.

    Args:
        pipeline (Pipeline): The OVITO pipeline object to render.
        imagefile (str): Path to save the rendered image (e.g., 'output.png').
        size (tuple): Image resolution as (width, height) in pixels.
        renderer (Renderer): OVITO renderer instance, e.g., TachyonRenderer.
        camera_dir (tuple): Camera direction vector.
        viewtype (Viewport.Type): Viewport projection type (e.g., orthographic or perspective).

    Example:
        >>> pipeline = import_file('POSCAR')
        >>> pipeline.modifiers.append(CommonNeighborAnalysisModifier())
        >>> my_render_image(
        ...     pipeline=pipeline,
        ...     imagefile='render.png',
        ...     size=(600, 2400),
        ...     renderer=TachyonRenderer(
        ...         shadows=True,
        ...         direct_light_intensity=1.2,
        ...         antialiasing_samples=64,
        ...         ambient_occlusion_samples=64,
        ...         max_ray_recursion=100,
        ...         depth_of_field=False
        ...     ),
        ...     camera_dir=(1, 0, 0)
        ... )
    """
    pipeline.add_to_scene()
    vp = Viewport()
    vp = Viewport(type = viewtype, camera_dir = camera_dir)
    vp.zoom_all(size=size)
    vp.render_image(filename=imagefile, 
                    size=size, 
                    renderer=renderer)
    pipeline.remove_from_scene()