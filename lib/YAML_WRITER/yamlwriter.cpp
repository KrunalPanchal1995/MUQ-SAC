#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <string>

namespace py = pybind11;

// Forward declaration to resolve circular dependency
YAML::Node convert_py_object_to_yaml(const py::handle& obj);

YAML::Node convert_py_dict_to_yaml(const py::dict& py_dict) {
    YAML::Node node;

    for (auto item : py_dict) {
        // Process each key-value pair in the dictionary
        std::string key = py::str(item.first); // Convert key to std::string

        // Convert value to YAML node
        YAML::Node yaml_value = convert_py_object_to_yaml(item.second);

        // Assign to YAML node with appropriate key
        node[key] = yaml_value;
    }

    return node;
}


// Function to convert a Python object to YAML node
YAML::Node convert_py_object_to_yaml(const py::handle& obj) {
    YAML::Node node;
    if (py::isinstance<py::bool_>(obj)) {
        node = obj.cast<bool>() ? "True" : "False";
    } else if (py::isinstance<py::int_>(obj)) {
        node = obj.cast<double>();
    } else if (py::isinstance<py::float_>(obj)) {
        node = obj.cast<double>();
    } else if (py::isinstance<py::str>(obj)) {
        node = obj.cast<std::string>();
    } else if (py::isinstance<py::tuple>(obj)) {
        const auto& tuple = obj.cast<py::tuple>();
        for (size_t i = 0; i < tuple.size(); ++i) {
            node.push_back(convert_py_object_to_yaml(tuple[i]));
        }
    } else if (py::isinstance<py::list>(obj)) {
        const auto& list = obj.cast<py::list>();
        for (size_t i = 0; i < list.size(); ++i) {
            node.push_back(convert_py_object_to_yaml(list[i]));
        }
    } else if (py::isinstance<py::dict>(obj)) {
        node = convert_py_dict_to_yaml(obj.cast<py::dict>());
    } else {
        throw std::runtime_error("Unsupported type in Python object");
    }
    return node;
}
// Function to dump a Python dictionary to a YAML file with custom path and filename
void dump_to_yaml(const std::string& path, const std::string& filename, const py::dict& data) {
    YAML::Node node = convert_py_dict_to_yaml(data);

    // Construct the full path including filename
    std::string filepath = path + "/" + filename;

    std::ofstream fout(filepath);
    fout << node;
}

// Pybind11 module definition
PYBIND11_MODULE(yamlwriter, m) {
    m.def("dump_to_yaml", &dump_to_yaml, "Dump dictionary to YAML file with custom path and filename",
          py::arg("path"), py::arg("filename"), py::arg("data"));
}

