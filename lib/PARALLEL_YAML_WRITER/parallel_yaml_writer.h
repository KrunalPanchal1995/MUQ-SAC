// parallel_yaml_writer.h
#ifndef PARALLEL_YAML_WRITER_H
#define PARALLEL_YAML_WRITER_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <iostream>
#include <fstream>
#include <thread>
#include <atomic>
#include <mutex>

namespace py = pybind11;
extern std::mutex mtx;

YAML::Node convert_py_dict_to_yaml(const py::dict &dict);
YAML::Node convert_py_object_to_yaml(const py::handle& obj);

void dump_to_yaml(const std::string &path, const std::string &filename, const py::dict &data);

void process_yaml_files(const std::vector<py::dict> &data, const std::vector<std::string> &path, const std::string &base_filename, int num_threads);

/*PYBIND11_MODULE(parallel_yaml_writer, m) {
    m.def("dump_to_yaml", &dump_to_yaml, "Dump a Python dictionary to a YAML file",
          py::arg("path"), py::arg("filename"), py::arg("data"));
    m.def("process_yaml_files", &process_yaml_files, "Process YAML files in parallel",
          py::arg("data"), py::arg("paths"), py::arg("base_filename"), py::arg("num_threads"));
}*/

#endif // PARALLEL_YAML_WRITER_H

