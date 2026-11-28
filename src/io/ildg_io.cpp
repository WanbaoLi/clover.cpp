// The MPI subarray layout follows Bridge++'s FieldIO_LIME_Parallel
// implementation by Tatsumi Aoyama (aoym).

#include "io/ildg_io.h"
#include "runtime/communicator_mpi.h"
#include "field/field_G.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check_mpi(int status, const char* operation)
{
  if (status == MPI_SUCCESS) return;
  char message[MPI_MAX_ERROR_STRING];
  int length = 0;
  MPI_Error_string(status, message, &length);
  throw std::runtime_error(std::string(operation) + ": " + std::string(message, length));
}

}  // namespace

std::uint64_t ildg_payload_offset(const std::string& path, const std::vector<int>& shape)
{
  std::ifstream file(path, std::ios::binary | std::ios::ate);
  if (!file) throw std::runtime_error("cannot open gauge: " + path);
  const auto end = file.tellg();
  if (end < 0) throw std::runtime_error("cannot determine gauge size: " + path);
  const auto bytes = static_cast<std::uint64_t>(end);

  constexpr std::uint64_t header_size = 144;
  std::uint64_t expected = 4 * 18 * sizeof(double);
  for (int extent : shape) {
    if (extent < 1 || expected > std::numeric_limits<std::uint64_t>::max() / extent)
      throw std::runtime_error("invalid ILDG lattice extent");
    expected *= extent;
  }

  std::uint64_t payload = 0;
  bool geometry = false;
  for (std::uint64_t at = 0; at < bytes;) {
    unsigned char header[header_size];
    file.seekg(at);
    file.read(reinterpret_cast<char*>(header), sizeof(header));
    if (!file || header[0] != 0x45 || header[1] != 0x67 || header[2] != 0x89 ||
        header[3] != 0xab || header[4] != 0 || header[5] != 1)
      throw std::runtime_error("invalid LIME record header");

    std::uint64_t length = 0;
    for (int i = 8; i < 16; ++i) length = (length << 8) | header[i];
    const auto content = at + header_size;
    if (length > bytes - content) throw std::runtime_error("truncated LIME record");

    const auto* record = header + 16;
    const auto* record_end = std::find(record, record + 128, 0);
    const std::string type(record, record_end);
    if (type == "ildg-binary-data") {
      if (payload || length != expected)
        throw std::runtime_error("expected one FP64 ILDG gauge payload of the declared lattice size");
      payload = content;
    } else if (type == "ildg-format" || type == "scidac-private-file-xml") {
      if (length > 4 * 1024 * 1024) throw std::runtime_error("unexpected XML record size");
      std::string xml(length, ' ');
      if (length && !file.read(&xml[0], length))
        throw std::runtime_error("truncated LIME XML record");

      std::smatch match;
      std::vector<int> found;
      const std::regex dims("<dims>\\s*([0-9 ]+(?:\\s+[0-9]+)*)\\s*</dims>");
      if (std::regex_search(xml, match, dims)) {
        std::istringstream scan(match[1].str());
        int extent;
        while (scan >> extent) found.push_back(extent);
      } else {
        for (const std::string tag : {"lx", "ly", "lz", "lt"}) {
          const std::regex axis("<" + tag + ">\\s*([0-9]+)\\s*</" + tag + ">");
          if (std::regex_search(xml, match, axis)) found.push_back(std::stoi(match[1]));
        }
      }
      if (!found.empty()) {
        if (found != shape) throw std::runtime_error("LIME geometry differs from --lattice");
        geometry = true;
      }
    }

    const auto padding = (8 - length % 8) % 8;
    if (padding > bytes - content - length) throw std::runtime_error("truncated LIME padding");
    at = content + length + padding;
  }
  if (!geometry || !payload)
    throw std::runtime_error("expected ILDG/SciDAC dimensions and ildg-binary-data");
  return payload;
}

void read_ildg_gauge(Field_G& field, const std::string& path, std::uint64_t payload)
{
  constexpr int dimensions = 4;
  auto comm = Communicator_impl::world();
  int global[] = {CommonParameters::Lx(), CommonParameters::Ly(), CommonParameters::Lz(), CommonParameters::Lt()};
  int local[] = {CommonParameters::Nx(), CommonParameters::Ny(), CommonParameters::Nz(), CommonParameters::Nt()};
  int starts[dimensions];
  for (int axis = 0; axis < dimensions; ++axis) starts[axis] = local[axis] * Communicator::ipe(axis);

  const int site_elements = field.nin() * field.nex();
  const int volume = field.nvol();
  std::vector<double> buffer(static_cast<std::size_t>(site_elements) * volume);
  MPI_Offset offset = static_cast<MPI_Offset>(payload);
  check_mpi(MPI_Bcast(&offset, 1, MPI_OFFSET, 0, comm), "broadcast ILDG payload offset");

  MPI_File file;
  check_mpi(MPI_File_open(comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file), "open ILDG file");
  MPI_Datatype site, subarray;
  check_mpi(MPI_Type_contiguous(sizeof(double) * site_elements, MPI_BYTE, &site), "create ILDG site type");
  check_mpi(MPI_Type_commit(&site), "commit ILDG site type");
  check_mpi(MPI_Type_create_subarray(dimensions, global, local, starts, MPI_ORDER_FORTRAN, site, &subarray),
            "create ILDG subarray");
  check_mpi(MPI_Type_commit(&subarray), "commit ILDG subarray");
  check_mpi(MPI_File_set_view(file, offset, site, subarray, "native", MPI_INFO_NULL), "set ILDG file view");
  MPI_Status status;
  check_mpi(MPI_File_read_all(file, buffer.data(), volume, site, &status), "read ILDG subarray");
  int read_sites = 0;
  check_mpi(MPI_Get_count(&status, site, &read_sites), "count ILDG sites read");
  check_mpi(MPI_File_close(&file), "close ILDG file");
  check_mpi(MPI_Type_free(&subarray), "release ILDG subarray");
  check_mpi(MPI_Type_free(&site), "release ILDG site type");
  if (read_sites != volume) throw std::runtime_error("short ILDG gauge read");

  const std::uint16_t one = 1;
  if (*reinterpret_cast<const unsigned char*>(&one) == 1) {
    auto* bytes = reinterpret_cast<unsigned char*>(buffer.data());
    const auto end = bytes + buffer.size() * sizeof(double);
    for (; bytes != end; bytes += sizeof(double)) std::reverse(bytes, bytes + sizeof(double));
  }

  auto value = buffer.begin();
  for (int index = 0; index < volume; ++index) {
    for (int direction = 0; direction < field.nex(); ++direction) {
      for (int component = 0; component < field.nin(); ++component)
        field.set(component, index, direction, *value++);
    }
  }
}
