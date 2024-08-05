from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools

__version__ = '0.0.1'

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to delay importing pybind11 until it is actually
    installed, so that the `get_include()` method can be invoked. """

    def __str__(self):
        import pybind11
        return pybind11.get_include()

ext_modules = [
    Extension(
        'shuffle',
        ['shuffle.cpp'],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
        ],
        language='c++'
    ),
]

setup(
    name='shuffle',
    version=__version__,
    author='Krunal Panchal',
    author_email='me18d005@iitm.ac.in',
    description='Shuffle arrays with multiple threads',
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.5.0'],
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)

