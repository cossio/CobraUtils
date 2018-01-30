"""
Input/output with the CHO_K1 model files.
"""

import os
import subprocess


def gunzip(path):
    """
    Check if file is gzipped, in which case extract it. 
    Returns original path.
    """
    assert os.path.exists(path) or os.path.exists(path + '.gz')
    if not os.path.exists(path):
        subprocess.call(['gunzip', '--keep', path + '.gz'])
    return path

