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

#ifndef BLOCK_READER_H
#define BLOCK_READER_H

#include <string>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <sstream>
#include <ios>

#include "basic_types.h"
#include "entry.h"

using namespace std;



class BlockReader {

public:

  BlockReader(const string &f, const Counter &m, const bool &u, const bool &que) {
    filename = f;
    threshold_cov = m;
    unweighted = u;
    unique = que;

    try {
      input.open(filename, ios::in);
    } catch(exception & e) { 
      cerr << "ERROR: failing opening the input file: " << filename << "\": " << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    if(!input.is_open()) {
      cerr << "ERROR: failing opening the input file: " << filename << endl;
      exit(EXIT_FAILURE);
    }

    already_got = false;
    end = false;
  }
  ~BlockReader() { }

  bool has_next() {
    if(unique) {
      return has_next_unique();
    } else {
      return has_next_nounique();
    }
  }

  Block get_block();

private:
  
  //Attributes
  string filename;
  Counter threshold_cov;
  bool unweighted;
  bool unique;

  ifstream input;

  bool end;
  bool already_got;
  vector<Fragment> fragment_block;
  Block block;
  unordered_set<Pointer> read_positions;
  Pointer max_position;
  Fragment last_fragment;
  vector<Fragment::const_iterator> fragment_pointers;

  //Private Methods
  bool has_next_unique();
  bool has_next_nounique();
  void extract_block();
  void string_to_fragment(const string &line, Fragment &read);
  void add_positions(const Fragment &read);
};

#endif

