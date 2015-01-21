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

#include "binomial.h"

std::vector<std::vector<unsigned int> > binom_coeff::btable;
std::vector<std::vector<unsigned int> > binom_coeff::ctable;

void
binom_coeff::initialize_binomial_coefficients(const unsigned int n,
                                              const unsigned int k) {
  // binomial coefficients
  btable.clear();
  btable.resize(n+1, std::vector<unsigned int>(n + 1, 0));
  for (unsigned int i = 1; i <= n; ++i) {
    for (unsigned int j = 0; j <= i; j++) {
      if (j == 0 || j == i) {
        btable[i][j] = 1;
      } else {
        btable[i][j] = btable[i - 1][j - 1] + btable[i - 1][j];
      }
    }
  }
  // cumulative binomial coefficients
  ctable.clear();
  ctable.resize(n+1, std::vector<unsigned int>(n + 1, 0));
  for (unsigned int i = 1; i <= n; i++) {
    for (unsigned int j = 0; j <= k; j++) {
      for(unsigned int x = 0; x <= j; x++) {
        ctable[i][j] += btable[i][x];
      }
    }
  }
}
