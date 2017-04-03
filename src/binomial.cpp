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
  for (unsigned int i = 0; i <= n; ++i) {
    for (unsigned int j = 0; j <= i; j++) {
      if (i == 0 || j == 0 || j == i) {
        btable[i][j] = 1;
      } else {
        btable[i][j] = btable[i - 1][j - 1] + btable[i - 1][j];
      }
    }
  }
  // cumulative binomial coefficients
  ctable.clear();
  ctable.resize(n+1, std::vector<unsigned int>(n + 1, 0));
  for (unsigned int i = 0; i <= n; i++) {
    for (unsigned int j = 0; j <= k; j++) {
      for(unsigned int x = 0; x <= j; x++) {
        ctable[i][j] += btable[i][x];
      }
    }
  }
}


unsigned int binom_coeff::indexof(std::bitset<MAX_COVERAGE> comb) {

  int k = 0;
  int c_k = 0;
  int temp = 0;
  int result = 0;

  while(comb.any()) {
    temp = ffsl(comb.to_ulong());
    c_k += temp;
    k++;
    result += binom_coeff::binomial_coefficient(c_k - 1, k);
    comb >>= (temp);
  }

  return result;
}


unsigned int binom_coeff::cumulative_indexof(std::bitset<MAX_COVERAGE> comb, const unsigned int n_elements) {

  unsigned int k = comb.count();
  unsigned int result = indexof(comb);

  for(unsigned int i = 0; i < k; i++) {
    result += binom_coeff::binomial_coefficient(n_elements, i);
  }

  return result;
}
