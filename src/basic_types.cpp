/**
 *
 *                              HapCol
 * Fast and Memory-efficient Haplotype Assembly From Gapless Long Reads
 *
 * Copyright (C) 2015  Yuri Pirola, Simone Zaccaria
 *
 * Distributed under the terms of the GNU General Public License (GPL)
 *
 * This file is part of HapCol.
 *
 * HapCol is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * HapCol is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HapCol.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#include "basic_types.h"

#include <limits>

const Cost::cost_t Cost::infinity_(std::numeric_limits<Cost::cost_t>::max());
const Cost Cost::INFTY(Cost::infinity_);


std::ostream& operator<<(std::ostream& out, const Cost& c) {
  if (c == Cost::INFTY)
    out << "INFINITY";
  else
    out << c.cost_;
  return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<bool>& v) {
  for (std::vector<bool>::const_iterator it= v.begin(); it != v.end(); ++it) {
    out << (*it? '1' : '0');
  }
  return out;
}


std::ostream& operator<<(std::ostream& out, const std::vector<char>& v) {
  for (std::vector<char>::const_iterator it= v.begin(); it != v.end(); ++it) {
    out << *it;
  }
  return out;
}



std::ostream& operator<<(std::ostream& out, const options_t& options) {
  const std::string SEP("\n");
  out
    << "Initialized? " << (options.options_initialized?"True":"False") << SEP
    << "Input filename: '" << options.input_filename << '\'' << SEP
    << "Haplotype filename: '" << options.haplotype_filename << '\'' << SEP
    << "Discard weights? " << (options.unweighted?"True":"False") << SEP
    << "Mask ambiguous positions? " << (options.no_xs?"False":"True") << SEP
    << "Error rate: " << options.error_rate << SEP
    << "Alpha: " << options.alpha;
  return out;
}

#include <boost/program_options.hpp>

options_t parse_arguments(int argc, char** argv) {
  namespace po = boost::program_options;

  options_t ret;
  po::options_description opts_desc("Program options");
  opts_desc.add_options()
    ("help,h", "produce (this) help message")
    ("input,i", po::value<std::string>(&ret.input_filename)->required(), "file containing the input reads (in WIF format)")
    ("haplotypes,o", po::value<std::string>(&ret.haplotype_filename)->required(), "file where the computed haplotypes will be written to")
    ("discard-weights,u", po::bool_switch(&ret.unweighted), "discard weights")
    ("no-ambiguous,x", po::bool_switch(&ret.no_xs), "do not mark ambiguous positions with Xs")
    ("error-rate,e", po::value<double>(&ret.error_rate)->default_value(ret.error_rate), "read error rate")
    ("alpha,a", po::value<double>(&ret.alpha)->default_value(ret.alpha), "significance (smaller is better)");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, opts_desc), vm);
    if (vm.count("help")) {
      std::cout << opts_desc << std::endl;
      return ret;
    }
    po::notify(vm);
    if ((ret.error_rate < 0.0) || (ret.error_rate > 1.0))  throw std::logic_error("error-rate must be a value between 0.0 and 1.0");
    if ((ret.alpha < 0.0) || (ret.alpha > 1.0))  throw std::logic_error("alpha must be a value between 0.0 and 1.0");
    ret.options_initialized= true;
  } catch (std::exception& e) {
    std::cout << "ERROR while parsing the program options: " << e.what() << std::endl;
    std::cout << opts_desc << std::endl;
  }

  return ret;
}
