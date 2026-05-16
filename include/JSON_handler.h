#pragma once
// Stub JSON_handler — always returns true.
// Replace with a real DCS JSON parser if run/lumi filtering is needed.
struct JSON_handler {
  bool isGood(unsigned int /*run*/, unsigned int /*lumi*/) const { return true; }
};
