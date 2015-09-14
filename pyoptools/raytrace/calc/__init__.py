'''Method collection to make calculation for optical systems
'''

from .calc import (chief_ray_search,
                   find_apperture,
                   find_ppp,
                   find_reference_sphere_radius,
                   get_optical_path_ep,
                   intersection,
                   nearest_points,
                   parallel_propagate,
                   parallel_propagate_ns,
                   paraxial_location,
                   pupil_location,
                   ray_paths)


__all__ = ["chief_ray_search",
           "find_apperture",
           "find_ppp",
           "find_reference_sphere_radius",
           "get_optical_path_ep",
           "intersection",
           "nearest_points",
           "parallel_propagate",
           "parallel_propagate_ns",
           "paraxial_location",
           "pupil_location",
           "ray_paths"]
