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

#ifndef _BASIC_TYPES_H_
#define _BASIC_TYPES_H_

#include <bitset>
#include <vector>
#include <iostream>
#include <limits>

#include "entry.h"

#define MAX_COVERAGE 32
#define MAX_CORRECTIONS 31

typedef unsigned int Counter;

#define MAX_COUNTER std::numeric_limits<Counter>::max()

typedef int Pointer;
typedef std::bitset<MAX_COVERAGE> BitColumn;
typedef std::vector<Entry> Column;


struct options_t {
  bool options_initialized;
  std::string input_filename;
  std::string haplotype_filename;
  bool unweighted;
  bool no_xs;
  double error_rate;
  double alpha;

  options_t()
  : options_initialized(false),
    input_filename(""),
    haplotype_filename(""),
    unweighted(true),
    no_xs(true),
    error_rate(0.05),
    alpha(0.01)
  {}

};

std::ostream& operator<<(std::ostream& out, const options_t& options);


struct constants_t
{

  BitColumn zeroes;
  BitColumn ones;

  constants_t() {
    ones.flip();
  };

};


struct Backtrace1
{
  Pointer jump;
  Counter index;

  Backtrace1()
    : jump(-1), index(0)
  {};
};



// A type for representing costs
class Cost {
public:
  typedef unsigned int cost_t;

private:
  static const cost_t infinity_;

  cost_t cost_;

  static bool is_addition_unsafe(const Cost& c1_, const Cost& c2_) {
    return (c1_.cost_ > (infinity_ - c2_.cost_));
  }

public:

  static const Cost INFTY;

  Cost(const cost_t cost= 0)
    :cost_(cost)
  {}

  Cost(const Cost& c)
    :cost_(c.cost_)
  {}

  Cost& operator=(const Cost& c) {
    cost_= c.cost_;
    return *this;
  }

  Cost& operator+=(const Cost& c) {
    if (is_addition_unsafe(*this, c))
      cost_ = infinity_;
    else
      cost_ += c.cost_;
    return *this;
  }

  Cost operator+(const Cost& c) const {
    if (is_addition_unsafe(*this, c))
      return INFTY;
    return Cost(cost_ + c.cost_);
  }

  bool operator<(const Cost& c) const {
    return cost_ < c.cost_;
  }

  bool operator<=(const Cost& c) const {
    return cost_ <= c.cost_;
  }

  bool operator>(const Cost& c) const {
    return cost_ > c.cost_;
  }

  bool operator>=(const Cost& c) const {
    return cost_ >= c.cost_;
  }

  bool operator==(const Cost& c) const {
    return cost_ == c.cost_;
  }

  friend std::ostream& operator<<(std::ostream& out, const Cost& c);
};

// Pretty-print costs
std::ostream& operator<<(std::ostream& out, const Cost& c);

// Pretty-print binary vectors
std::ostream& operator<<(std::ostream& out, const std::vector<bool>& v);

// Pretty-print char vectors
std::ostream& operator<<(std::ostream& out, const std::vector<char>& v);


options_t parse_arguments(int argc, char** argv);


#endif // _BASIC_TYPES_H_
