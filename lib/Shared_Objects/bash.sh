#!/bin/bash
python3.9 setup.py build_ext --inplace
python3.9 rename_so.py ./ ./
