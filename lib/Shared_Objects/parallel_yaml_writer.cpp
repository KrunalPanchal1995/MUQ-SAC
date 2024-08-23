#include "parallel_yaml_writer.h"
std::atomic<int> progress(0);
std::mutex mtx;
auto start_time = std::chrono::high_resolution_clock::now();


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
        // Process each key-value pair in the dictionary
        std::string key = py::str(item.first); // Convert key to std::string
	//std::cout << "Key: " << key;
        // Convert value to YAML node
        //std::cout << " Key: " << key<<" Value: "<<item.second<<std::endl;
        YAML::Node yaml_value = convert_py_object_to_yaml(item.second);
	
        // Assign to YAML node with appropriate key
        node[key] = yaml_value;
        //std::cout << " Key: " << key<<" Value: "<<yaml_value<<std::endl;
    }

    return node;
}

YAML::Node convert_py_object_to_yaml(const py::handle& obj) {
    YAML::Node node;
    //py::handle bool_type = py::module::import("builtins").attr("bool");
    //type t = obj.get_type();
    //std::cout << "Value type: " << t <<std::endl;
    if (py::isinstance<py::bool_>(obj)) {  // Handle boolean values
        //std::cout << "Value: bool" << obj;
        node = obj.cast<bool>() ? "True" : "False"; // Convert bool to string "true" or "false"
    } else if (py::isinstance<py::int_>(obj)) {
        // Convert int to double to preserve high precision numbers
        //std::cout << "Value: int" << obj;
        node = obj.cast<double>();
    } else if (py::isinstance<py::float_>(obj)) {
        // Convert float to double to preserve high precision numbers
        node = obj.cast<double>();
    } else if (py::isinstance<py::str>(obj)) {
        //std::cout << "Value: string" << obj;
        node = obj.cast<std::string>();
    } else if (py::isinstance<py::tuple>(obj)) {
        //std::cout << "Value: tuple" << obj;
        const auto& tuple = obj.cast<py::tuple>();
        for (size_t i = 0; i < tuple.size(); ++i) {
            node.push_back(convert_py_object_to_yaml(tuple[i]));
        }
    } else if (py::isinstance<py::list>(obj)) {
        //std::cout << "Value: list" << obj;
        const auto& list = obj.cast<py::list>();
        for (size_t i = 0; i < list.size(); ++i) {
            node.push_back(convert_py_object_to_yaml(list[i]));
        }
    } else if (py::isinstance<py::dict>(obj)) {
        //std::cout << "Value: dict" << obj;
        node = convert_py_dict_to_yaml(obj.cast<py::dict>());
    } else {
        throw std::runtime_error("Unsupported type in Python object");
    }

    return node;
}

void dump_to_yaml(const std::string& path, const std::string& filename, const py::dict& data) {
    std::lock_guard<std::mutex> lock(mtx);
    YAML::Node node = convert_py_dict_to_yaml(data);

    // Construct the full path including filename
    std::string filepath = path + "/" + filename;

    std::ofstream fout(filepath);
    fout << node;
    fout.close();
}


void process_chunk(const std::vector<py::dict>& data, const std::vector<std::string>& paths, const std::string& base_filename, int start, int end) {
    for (int i = start; i < end; ++i) {
        std::string filename = base_filename + ".yaml";
        dump_to_yaml(paths[i], filename, data[i]);

        // Update progress
        {
            std::lock_guard<std::mutex> lock(mtx);
            ++progress;
            float percent = static_cast<float>(progress) / data.size() * 100;
            std::cout << "\rProgress: " << percent << "% " << std::flush;
        }
    }
}

void process_yaml_files(const std::vector<py::dict>& data, const std::vector<std::string>& paths, const std::string& base_filename, int num_threads) {
    std::vector<std::thread> threads;
    int chunk_size = (data.size() + num_threads - 1) / num_threads;

    // Start threads
    for (int i = 0; i < num_threads; ++i) {
        int start = i * chunk_size;
        int end = std::min(start + chunk_size, static_cast<int>(data.size()));
        threads.emplace_back(process_chunk, std::ref(data), std::ref(paths), std::ref(base_filename), start, end);
    }

    // Join threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Clear progress bar after completion
    std::cout << "\rProgress: 100%   " << std::endl;
    // Calculate and print total time taken in hours, minutes, and seconds
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
          py::arg("data"), py::arg("paths"), py::arg("base_filename"), py::arg("num_threads"));
}

