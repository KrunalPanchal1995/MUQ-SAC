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
#include <algorithm>

namespace py = pybind11;

std::atomic<int> progress(0);
std::mutex mtx;
auto start_time = std::chrono::high_resolution_clock::now();

YAML::Node convert_py_dict_to_yaml(const py::dict& py_dict);
YAML::Node convert_py_object_to_yaml(const py::handle& obj);

namespace YAML {

// Convert pybind11::handle to YAML::Node
template<>
struct convert<py::handle> {
    static Node encode(const py::handle& rhs) {
        return convert_py_object_to_yaml(rhs);
    }

    static bool decode(const Node& node, py::handle& rhs) {
        // Optional: Implement decoding if needed
        return false;  // Example implementation; customize as needed
    }
};

}  // namespace YAML

YAML::Node convert_py_dict_to_yaml(const py::dict& py_dict) {
    YAML::Node node;
    for (auto item : py_dict) {
        std::string key = py::str(item.first);  // Convert key to std::string
        YAML::Node yaml_value = convert_py_object_to_yaml(item.second);
        node[key] = yaml_value;
    }
    return node;
}

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

void dump_to_yaml(const std::string& path, const std::string& filename, const py::dict& data) {
    std::lock_guard<std::mutex> lock(mtx);
    YAML::Node node = convert_py_dict_to_yaml(data);
    std::string filepath = path + "/" + filename;
    std::ofstream fout(filepath);
    if (!fout.is_open()) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }
    fout << node;
    fout.close();
}

void process_chunk(const std::vector<py::dict>& data, const std::vector<std::string>& paths, const std::string& base_filename, const std::vector<std::string>& indices, int start, int end) {
    for (int i = start; i < end; ++i) {
        std::string filename = base_filename + indices[i] + ".yaml";
        dump_to_yaml(paths[i], filename, data[i]);
        {
            std::lock_guard<std::mutex> lock(mtx);
            ++progress;
            float percent = static_cast<float>(progress) / data.size() * 100;
            std::cout << "\rProgress: " << percent << "% " << std::flush;
        }
    }
}

void process_yaml_files(const std::vector<py::dict>& data, const std::vector<std::string>& paths, const std::string& base_filename, const std::vector<std::string>& indices, int num_threads) {
    std::vector<std::thread> threads;
    using size_type = std::vector<py::dict>::size_type;
    size_type data_size = data.size();
    int chunk_size = 10;//(data_size + num_threads - 1) / num_threads;
	
	
    for (const auto& key: data){
        if (threads.size() >= num_threads) {
            for (auto& t : threads) {
                t.join();
            }
            threads.clear();
        }
        
        threads.emplace_back(process_chunk, std::ref(data), std::ref(paths), std::ref(base_filename), std::ref(indices), chunk_start, chunk_end);
    }
    
    
    //for (int start = 0; start < static_cast<int>(data_size); start += chunk_size) {
    //    for (int i = 0; i < num_threads && start + i * chunk_size < static_cast<int>(data_size); ++i) {
   //         int chunk_start = start + i * chunk_size;
   //         int chunk_end = std::min(chunk_start + chunk_size, static_cast<int>(data_size));
   //         threads.emplace_back(process_chunk, std::ref(data), std::ref(paths), std::ref(base_filename), std::ref(indices), chunk_start, chunk_end);
   //     }

   //     for (auto& thread : threads) {
   //         thread.join();
   //     }
   //     threads.clear();
   // }

    std::cout << "\rProgress: 100%   " << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    int total_seconds = static_cast<int>(elapsed_seconds.count());
    int hours = total_seconds / 3600;
    int minutes = (total_seconds % 3600) / 60;
    int seconds = total_seconds % 60;

    std::cout << "\nTime taken: " << hours << "h ";
    if (minutes < 10) std::cout << "0";
    std::cout << minutes << "m ";
    if (seconds < 10) std::cout << "0";
    std::cout << seconds << "s" << std::endl;
}

PYBIND11_MODULE(parallel_yaml_writer, m) {
    m.def("dump_to_yaml", &dump_to_yaml, "Dump a Python dictionary to a YAML file",
          py::arg("path"), py::arg("filename"), py::arg("data"));
    m.def("process_yaml_files", &process_yaml_files, "Process YAML files in parallel",
          py::arg("data"), py::arg("paths"), py::arg("base_filename"), py::arg("indices"), py::arg("num_threads"));
}

