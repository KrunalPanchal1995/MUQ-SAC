from setuptools import setup, Extension, Command
import platform
import shutil
import os,re
import sys
import subprocess
import stat
import glob
# Ensure pybind11 is installed
try:
    import pybind11
except ImportError:
    import subprocess
    import sys
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pybind11'])
    import pybind11

# Get pybind11 include directory
pybind11_include_dir = pybind11.get_include()

# Specify additional include directories if needed (e.g., where your .h files are located)
include_dirs = [pybind11_include_dir]
home_path = os.getcwd()
# List of C++ files and header files
cpp_files = ['parallel_yaml_writer.cpp', 'shuffle.cpp','yamlwriter.cpp']
header_files = ['parallel_yaml_writer.h']  # Add your header files here

# Combine both C++ and header files for permissions adjustment
all_files = cpp_files + header_files

# Adjust file permissions
for file in all_files:
    if platform.system() in ['Linux', 'Darwin']:  # macOS and Ubuntu (Linux)
        os.chmod(file, os.stat(file).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    elif platform.system() == 'Windows':
        # On Windows, we generally don't need to change file permissions
        pass

# Function to generate module names and output filenames from cpp file names
def generate_module_name(filename):
    base_name = os.path.splitext(filename)[0]  # Remove .cpp extension
    return re.sub(r'[^a-zA-Z0-9_]', '_', base_name)  # Replace non-alphanumeric characters with underscores

# Function to generate the output filename for .so files
def generate_so_filename(filename):
    module_name = generate_module_name(filename)
    if platform.system() == 'Windows':
        return f"{module_name}.pyd"  # Windows uses .pyd for Python extensions
    else:
        return f"lib{module_name}.so"  # Linux/macOS uses .so for shared libraries

# Define the extensions
extensions = [
    Extension(
        name=generate_module_name(cpp_file),
        sources=[cpp_file],
        include_dirs=include_dirs,  # Include directories for headers
        extra_compile_args=['-std=c++11'] if platform.system() != 'Windows' else [],
        extra_link_args=[] if platform.system() != 'Windows' else []
    )
    for cpp_file in cpp_files
]

# Requirements that need to be installed
#requirements = [
#    'pybind11>=2.2',
#    'numpy>=1.18.0',
    # Add other Python package dependencies here
#]

## OR

# Read requirements from requirements.txt
with open('requirements.txt') as f:
    requirements = f.read().splitlines()
# Setup configuration

class RenameBuildExtCommand(Command):
    """Custom command to rename .so files after building."""
    description = 'Rename shared object files to remove platform-specific suffixes'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        # Run the standard build_ext command
        self.run_command('build_ext')
        
        build_lib = os.path.abspath(self.distribution.get_command_obj('build_ext').build_lib)
        for cpp_file in cpp_files:
            module_name = generate_module_name(cpp_file)
            old_pattern = os.path.join(build_lib, f"{module_name}.cpython-*.so")
            new_filename = os.path.join(build_lib, f"{module_name}.so")
            
            for file in glob.glob(old_pattern):
                if file != new_filename:  # Avoid renaming if already correct
                    shutil.move(file, new_filename)
                    print(f"Renamed {file} to {new_filename}")


setup(
    name='cpp_extension_module',
    version='1.0',
    description='A package with C++ extensions',
    ext_modules=extensions,
    install_requires=requirements,  # Directly specified requirements
    setup_requires=['setuptools', 'pybind11>=2.2'],
    cmdclass={
        'rename_build_ext': RenameBuildExtCommand,
    }

)

