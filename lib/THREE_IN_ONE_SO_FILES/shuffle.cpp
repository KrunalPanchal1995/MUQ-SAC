#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include <random>
#include <thread>
#include <mutex>
#include <algorithm>

namespace py = pybind11;

std::mutex g_mutex;

void shuffle_array(py::array_t<double> array, int num_shuffles) {
    auto buf = array.request();
    double* ptr = (double*)buf.ptr;
    std::vector<double> data(ptr, ptr + buf.size);

    std::random_device rd;
    std::mt19937 g(rd());

    for (int i = 0; i < num_shuffles; ++i) {
        std::shuffle(data.begin(), data.end(), g);
    }

    std::lock_guard<std::mutex> guard(g_mutex);
    std::copy(data.begin(), data.end(), ptr);
}

void shuffle_arrays(py::dict V_, int num_threads, std::vector<std::string> unsrt, int sim) {
    int num_shuffles = int(sim * 0.8);
    std::vector<std::thread> threads;

    for (const auto& key : unsrt) {
        if (threads.size() >= num_threads) {
            for (auto& t : threads) {
                t.join();
            }
            threads.clear();
        }

        py::array_t<double> array = V_[key.c_str()].cast<py::array_t<double>>();
        threads.emplace_back(std::thread(shuffle_array, array, num_shuffles));
    }

    for (auto& t : threads) {
        t.join();
    }
}

PYBIND11_MODULE(shuffle, m) {
    m.def("shuffle_arrays", &shuffle_arrays, "Shuffle arrays with multiple threads",
          py::arg("V_"), py::arg("num_threads"), py::arg("unsrt"), py::arg("sim"));
}

