from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools

__version__ = '0.0.1'

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to delay importing pybind11 until it is actually
    installed, so that the `get_include()` method can be invoked. """
    
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

ext_modules = [
    Extension(
        'parallel_yaml_writer',
        ['parallel_yaml_writer.cpp'],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            # Path to yaml-cpp headers
            '/usr/local/include'  # Adjust this path based on your yaml-cpp installation
        ],
        libraries=['yaml-cpp'],
        library_dirs=[
            '/usr/local/lib'  # Adjust this path based on your yaml-cpp installation
        ],
        language='c++'
    ),
]

setup(
    name='parallel_yaml_writer',
    version=__version__,
    author='Krunal Panchal',
    author_email='me18d005@iitm.ac.in',
    description='Parallel YAML writer with multiple threads',
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.5.0'],
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)

