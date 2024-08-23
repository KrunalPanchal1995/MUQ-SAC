import os
import re
import glob
import shutil
import sys

def generate_module_name(filename):
    base_name = os.path.splitext(filename)[0]  # Remove extension
    return re.sub(r'[^a-zA-Z0-9_]', '_', base_name)  # Replace non-alphanumeric characters with underscores

def rename_so_files(build_lib,setup_lib):
    # Rename .so files
    for file in glob.glob(os.path.join(build_lib, '*.so')):
        base_name = os.path.basename(file)
        file_name = base_name.split('.')[0]
        #raise AssertionError("Stop")
        new_name = generate_module_name(file_name) + '.so'
        new_path = os.path.join(setup_lib, new_name)
        if file != new_path:  # Avoid renaming if already correct
            shutil.move(file, new_path)
            print(f"Renamed {file} to {new_path}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python rename_so_files.py <build_lib> <final_location_dir>")
        sys.exit(1)
    
    build_lib = sys.argv[1]
    setup_lib = sys.argv[2]
    rename_so_files(build_lib,setup_lib)

