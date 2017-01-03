#include "balanced_combinations.h"
#include <cmath>

using namespace std;


BalancedCombinations::BalancedCombinations() : generator() {}


void BalancedCombinations::initialize(const Counter n, const Counter k,
				      BitColumn col, const double r) {

  n_ = n;
  k_ = k;
  col_ = col;
  r_ = r;

  c_ = (Counter)ceil(n_ * r_);

  // pi_0 and pi_1
  p.push_back(n_ - col_.count());
  p.push_back(col_.count());

  // ..
  #ifdef BALANCED_COMBINATIONS_DEBUG
  cout << "n, k, c : " << n_ << ", " << k_ << ", " << c_ << endl;
  cout << "p0, p1 : " << p[0] << ", " << p[1] << " ..." << endl;
  #endif
  // ..

  // build mapping for composing combinations and initialize arrays
  build_mapping();
  initialize_arrays();

  // ..
  #ifdef BALANCED_COMBINATIONS_DEBUG
  all_combinations();
  #endif
  // ..

  // initialize the counters
  t_ = 0;
  i_ = 0;
  j_ = 0;
  ii_ = 0;
  jj_ = 0;

  // ..
  #ifdef BALANCED_COMBINATIONS_DEBUG
  cout << endl << "<<-- BalancedCombinations::try_next() -->>" << endl;
  #endif
  // ..

  has_next_ = true;
  s_ = true; // prime the try loop
  try_next();
}


bool BalancedCombinations::has_next() {

  return has_next_;
}


void BalancedCombinations::next() {

  make_current();
  s_ = false;
  try_next();
}


void BalancedCombinations::get_combination(BitColumn & result) {

  result.reset();
  result |= current_;
}


// auxiliary (private) functions
/**********************************************************************/


void BalancedCombinations::build_mapping() {

  map.resize(2);
  for(i_ = 0; i_ < n_; ++i_) {

    if(col_.test(i_))
      map[1].push_back(i_);
    else
      map[0].push_back(i_);
  }

  // ..
  #ifdef BALANCED_COMBINATIONS_DEBUG
  cout << endl << "<<-- BalancedCombinations::build_mapping() -->>" << endl;
  cout << "map_0 : " << vector_to_string(map[0]) << endl;
  cout << "map_1 : " << vector_to_string(map[1]) << endl;
  #endif
  // ..
}


void BalancedCombinations::initialize_arrays() {

  // c[0][.]
  a.clear();
  a.resize(p[0]+1);
  c.push_back(a);

  // c[1][.]
  a.clear();
  a.resize(p[1]+1);
  c.push_back(a);
}


void BalancedCombinations::retrieve_c0() {

  if(c[0][i_].empty()) {

    generator.initialize(p[0], i_);
    while(generator.has_next()) {

      generator.next(); // should always be at least the empty comb
      generator.get_combination(comb);
      c[0][i_].push_back(comb);

      // ..
      #ifdef BALANCED_COMBINATIONS_DEBUG
      cout << "add to c0 : " << comb << endl;
      #endif
      // ..

    }
  }
}


void BalancedCombinations::retrieve_c1() {

  if(c[1][j_].empty()) {

    generator.initialize(p[1], j_);
    while(generator.has_next()) {

      generator.next(); // should always be at least the empty comb
      generator.get_combination(comb);
      c[1][j_].push_back(comb);

      // ..
      #ifdef BALANCED_COMBINATIONS_DEBUG
      cout << "add to c1 : " << comb << endl;
      #endif
      // ..

    }
  }
}


void BalancedCombinations::make_current() {

  current_.reset();

  // fill c0
  for(i = 0; i < p[0]; ++i)
    if(c[0][i_][ii_].test(i))
      current_.set(map[0][i]);

  // fill c1
  for(j = 0; j < p[1]; ++j)
    if(c[1][j_][jj_].test(j))
      current_.set(map[1][j]);
}


void BalancedCombinations::try_next() {

  // loop with switch, advancing exactly once each call to function
  while(t_ <= k_) {
    while(i_ <= min(p[0], t_)) {
      j_ = t_ - i_;

      // check if j_ is feasible
      if(j_ <= p[1]) {

	// check balance threshold
	if((p[0]-i_ + min(p[1],t_-i_) >= c_) and (p[1]-j_ + min(p[0],t_-j_) >= c_)) {

	  retrieve_c0(); // c[0][i_]
	  while(ii_ < c[0][i_].size()) {

	    retrieve_c1(); // c[1][j_]
	    while(jj_ < c[1][j_].size()) {

	      // ..
              #ifdef BALANCED_COMBINATIONS_DEBUG
	      cout << "\t\t\t t, i, j, ii, jj : " << t_ << ", " << i_ << ", " << j_;
	      cout << ", " << ii_ << ", " << jj_ << endl;
	      cout << ". s is " << (s_ ? "true" : "false" ) << endl;
              #endif
	      // ..

	      // at this point, jj_,ii_,j_,i_,t_ is a valid configuration
	      if(s_)
		return;

	      s_ = true;
	      ++jj_;
	    }
	    jj_ = 0;
	    ++ii_;
	  }
	  ii_ = 0;
	}
      }
      ++i_;
    }
    i_ = 0;
    ++t_;
  }

  has_next_ = false; // the end
}


// for debugging ...
/**********************************************************************/
#ifdef BALANCED_COMBINATIONS_DEBUG


string BalancedCombinations::column_to_string(const BitColumn & col,
					      const Counter & len) {

  string str = col.to_string();
  reverse(str.begin(), str.end());

  return str.substr(0, len);
}


string BalancedCombinations::vector_to_string(const vector<Counter> & v) {

  string str = "";
  for(i = 0; i < v.size(); ++i)
    str.append(to_string(v[i]) + " ");

  return str;
}


void BalancedCombinations::all_combinations(bool verbose) {

  cout << endl << "<<-- BalancedCombinations::all_combinations() -->>" << endl;
  for(t_ = 0; t_ <= k_; ++t_) {
    cout << endl << "t : " << t_ << " -------------------------------" << endl;

    for(i_ = 0; i_ <= min(p[0], t_); ++i_) {
      j_ = t_ - i_;

      // check if j_ is feasible
      if(j_ > p[1])
	continue;

      cout << "i, j : " << i_ << ", " << j_ << " ..." << endl;

      // check balance threshold
      if((p[0] - i_ + min(p[1], t_ - i_)) < c_)
	continue;
      if((p[1] - j_ + min(p[0], t_ - j_)) < c_)
	continue;

      retrieve_c0(); // c[0][i_]
      for(ii_ = 0; ii_ < c[0][i_].size(); ++ii_) {

	retrieve_c1(); // c[1][j_]
	for(jj_ = 0; jj_ < c[1][j_].size(); ++jj_) {

	  if(verbose) {
	    cout << endl << c[0][i_][ii_];
	    cout << " + " << c[1][j_][jj_] << endl;
	  }
	  else {
	    cout << "  " << column_to_string(c[0][i_][ii_], p[0]) << " + ";
	    cout << column_to_string(c[1][j_][jj_], p[1]);
	  }

	  make_current(); // current_ should hold the current value

	  cout << " -> ";
	  if(verbose)
	    cout << current_;
	  else
	    cout << column_to_string(current_, n_);
	  cout << endl;
	}
      }
    }
  }
}

#endif // BALANCED_COMBINATIONS_DEBUG
