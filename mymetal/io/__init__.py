"""
mymetal.io

This subpackage provides input and output utilities for interacting with VASP files 
and post-processing functions used in material simulations.

It includes functions for reading and writing VASP input/output files, as well 
as creating and writing content to files for post-processing tasks.

Modules:
    - vasp: Module for reading and writing VASP files.
    - post: Module for post-processing tasks.
"""

from mymetal.io.vasp import my_read_vasp, my_write_vasp, read_vasp_list
from mymetal.io.post.construct import create_content
from mymetal.io.post.write import write_content_to_file

__all__ = ['my_read_vasp', 'my_write_vasp', 'read_vasp_list',
            'create_content', 
            'write_content_to_file'
]
