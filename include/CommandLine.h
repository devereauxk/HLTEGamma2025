#pragma once
#include <string>
#include <map>
#include <cstdlib>

class CommandLine {
  std::map<std::string, std::string> args_;
public:
  CommandLine(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
      std::string a = argv[i];
      if (a.size() > 2 && a[0] == '-' && a[1] == '-' && i + 1 < argc)
        args_[a.substr(2)] = argv[++i];
    }
  }
  std::string Get(const std::string &k, const std::string &d = "") const {
    auto it = args_.find(k);
    return it != args_.end() ? it->second : d;
  }
  int GetInt(const std::string &k, int d = 0) const {
    auto it = args_.find(k);
    return it != args_.end() ? std::stoi(it->second) : d;
  }
  double GetDouble(const std::string &k, double d = 0) const {
    auto it = args_.find(k);
    return it != args_.end() ? std::stod(it->second) : d;
  }
  bool GetBool(const std::string &k, bool d = false) const {
    auto it = args_.find(k);
    if (it == args_.end()) return d;
    return it->second == "true" || it->second == "1";
  }
};
