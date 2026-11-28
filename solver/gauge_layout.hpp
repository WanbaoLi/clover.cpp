#pragma once
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Bridge++'s parallel reader uses the caller's geometry. Check the LIME header
// before allocating fields so a different shape with the same volume is rejected.
inline void validate_lime(const std::string& path,const std::vector<int>& shape) {
  std::ifstream f(path.c_str(),std::ios::binary|std::ios::ate);
  if (!f) throw std::runtime_error("cannot open gauge: "+path);
  const std::uint64_t bytes=static_cast<std::uint64_t>(f.tellg());
  std::uint64_t at=0,expected=4*18*8;
  for (int n:shape) expected*=n;
  bool binary=false,geometry=false;
  while (at<bytes) {
    unsigned char h[144];f.seekg(at);f.read(reinterpret_cast<char*>(h),sizeof(h));
    if (!f || h[0]!=0x45 || h[1]!=0x67 || h[2]!=0x89 || h[3]!=0xab || h[4]!=0 || h[5]!=1)
      throw std::runtime_error("invalid LIME record header");
    std::uint64_t len=0;
    for (int i=8;i<16;++i) len=(len<<8)|h[i];
    if (len>bytes-at-144) throw std::runtime_error("truncated LIME record");
    const std::string type(reinterpret_cast<const char*>(h+16),std::find(h+16,h+144,0)-(h+16));
    if (type=="ildg-binary-data") {
      if (binary || len!=expected) throw std::runtime_error("expected one FP64 ILDG gauge payload of the declared lattice size");
      binary=true;
    } else if (type=="ildg-format" || type=="scidac-private-file-xml") {
      if (len>4*1024*1024) throw std::runtime_error("unexpected XML record size");
      std::string xml(len,' ');if (len) f.read(&xml[0],len);
      std::smatch match;std::vector<int> found;
      if (std::regex_search(xml,match,std::regex("<dims>\\s*([0-9 ]+(?:\\s+[0-9]+)*)\\s*</dims>"))) {
        std::istringstream scan(match[1].str());int n;
        while (scan>>n) found.push_back(n);
      } else {
        for (const std::string tag:{"lx","ly","lz","lt"}) {
          if (std::regex_search(xml,match,std::regex("<"+tag+">\\s*([0-9]+)\\s*</"+tag+">"))) found.push_back(std::stoi(match[1]));
        }
      }
      if (!found.empty()) {
        if (found!=shape) throw std::runtime_error("LIME geometry differs from --lattice");
        geometry=true;
      }
    }
    const auto pad=(8-len%8)%8;
    if (pad>bytes-at-144-len) throw std::runtime_error("truncated LIME padding");
    at+=144+len+pad;
  }
  if (!geometry || !binary) throw std::runtime_error("expected ILDG/SciDAC dimensions and ildg-binary-data");
}
