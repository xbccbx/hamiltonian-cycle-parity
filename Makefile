CXX = g++
CXXFLAGS ?= -O3 -std=c++17
OPENMP_FLAGS ?= -fopenmp
BUILD_DIR := build

.PHONY: all check
all: $(BUILD_DIR)/bench_poly $(BUILD_DIR)/bench_dp $(BUILD_DIR)/extended

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/bench_poly: artifacts/theory/code/bench_poly.cpp artifacts/theory/code/baseline.cpp artifacts/theory/code/input_families.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

$(BUILD_DIR)/bench_dp: artifacts/theory/code/baseline.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

$(BUILD_DIR)/extended: artifacts/practical32/code/extended.cpp artifacts/practical32/code/baseline.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

check: all
	./$(BUILD_DIR)/bench_poly verify
	./$(BUILD_DIR)/extended verify
	./$(BUILD_DIR)/bench_dp selftest
