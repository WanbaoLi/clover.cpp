#pragma once

#include <cstdint>
#include <string>
#include <vector>

class Field_G;

// Inspect the record once on MPI rank zero, before allocating gauge fields.
std::uint64_t ildg_payload_offset(const std::string& path, const std::vector<int>& shape);
// All ranks call this function. Only rank zero supplies the inspected offset.
void read_ildg_gauge(Field_G& field, const std::string& path, std::uint64_t payload);
