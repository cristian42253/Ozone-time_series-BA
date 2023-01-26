#!/usr/bin/env python3

""" 
"""
import os
from baycntsa import processing # settings, drawler, functions, processing, main

__version__ = '0.0.21'
__license__ = 'GNU'
__author__ = 'Cristian E. Garcia-Bermudez. <cristian.garcia@postgrado.uv.cl>'

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

__all__ = [
    "processing"
    # "settings", 
    # "drawler", 
    # "functions", 
    # "processing", 
    # "main"
]