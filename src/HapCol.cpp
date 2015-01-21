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

#include <stdlib.h>
#include <string.h>
#include <limits>
#include <math.h>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <stdexcept>


#include "basic_types.h"
#include "binomial.h"
#include "combinations.h"
#include "columnreader.h"

#ifdef LOAD_REVISION
#include "revision.h"
#endif

// Log messages with DEBUG priority and higher
#define LOG_MSG
#define LOG_THRESHOLD LOG_LEVEL_INFO
// Include log facilities. It should the last include!!
#include "log.h"




using namespace std;


static inline
Pointer next(const Pointer &indexer_pointer, const int &total_size, const int &shift)
{
  return (indexer_pointer + shift) % total_size;
}

static inline
Pointer prev(const Pointer &indexer_pointer, const int &total_size, const int &shift)
{
  return (indexer_pointer + total_size - shift) % total_size;
}

static inline
bool check_end(ColumnReader &column_reader, const vector<Column> &input, const Pointer &pointer)
{
  return !column_reader.has_next() && (input[pointer][0].get_read_id() == -1);
}

static inline
void complement_mask(BitColumn &mask, const Counter &length, const constants_t &constants)
{
  mask ^= (constants.ones<<length).flip();
}

static inline
string column_to_string(const BitColumn &mask, const unsigned int &len) {
  string str = mask.to_string();
  reverse(str.begin(), str.end());
  return str.substr(0, len);
}


template <typename T>
static inline
void replace_if_less(T& a, const T& b) {
  if(a > b) { a = b;}
}


void computeInputParams(Counter &num_cols, Counter &MAX_COV, Counter &MAX_L,
                        Counter &MAX_K, vector<Counter> &sum_successive_L,
                        Counter &coverage, const string &inputfilename, const bool &unweighted,
                        vector<vector<Counter> > &scheme_backtrace);
void intersect(const Column &colQ, const Column &colJ, const int &q,
               vector<vector<Pointer> > &forw_indexer, vector<vector<Pointer> > &back_indexer);
void represent_column(const Column &column, BitColumn &result, Counter &cov);
void cut(BitColumn &to_be_cut, const vector<Pointer> &indexer, Counter &active_pj);
void extract_common_mask(const Column &column_q, const Pointer &q_pointer,
                         const Column &column_j, const BitColumn &mask_colj,
                         const vector<vector<Pointer> > &back_indexer,
                         const vector<vector<Pointer> > &forw_indexer,
                         BitColumn &mask_qj, Counter &active_qj);
int compute_active_common(const Column &colJ, const Column &colQ);
void insert_col_and_update(vector<Column> &input, vector<Counter> &k_j, vector <Counter> &homo_cost,
                           vector<Cost> &homo_weight, const Pointer &pointer, const Column &column, const bool &unweighted,
                           vector<bool> &kind_homozygous, const Counter &step);
void compute_weight_mask(const BitColumn &mask, const Column &column, Cost &weight_mask);
void reconstruct_haplotypes(const vector<vector<vector<Backtrace1> > > &backtrace_table1,
                            const vector<vector<vector<bool> > > &backtrace_table2_haplotypes,
                            const vector<vector<vector<bool> > > &backtrace_table2_new_block,
                            const vector<bool> &is_homozygous,
                            const vector<bool> &homo_haplotypes,
                            const vector<Backtrace1> &best_heterozygous1,
                            const vector<bool> &best_heterozygous2_haplotypes,
                            const vector<bool> &best_heterozygous2_new_block,
                            vector<bool> &haplotype1, vector<bool> &haplotype2);
Counter computeK(const Counter &cov, const double &alpha = 0.0, const double &error_rate = 0.0);
void add_xs(const vector<bool> &haplo1, const vector<bool> &haplo2, 
            vector<char> &haplo1_out, vector<char> &haplo2_out,
            const string &input_filename, const unsigned int coverage, const bool unweighted);
int map_fragment(const vector<bool> &read, const vector<unsigned int> &weights, const unsigned int offset,
                 const vector<bool> &haplo1_in, const vector<bool> &haplo2_in, unsigned int &total_errors);
void make_haplo(const vector<bool> &haplo, const vector<vector<bool> > &mapping_haplo, vector<char> &haplo_out);
void count_alleles(const vector<bool> &col, vector<int> &counter);






int main(int argc, char** argv)
{
#if defined(VCS_DATE) && defined(VCS_SHORT_HASH) && defined(VCS_WC_MODIFIED)
  INFO("HapCol (" VCS_BRANCH "@" VCS_SHORT_HASH << (VCS_WC_MODIFIED ? "-dirty" : "-clean") << ")");
#else
  INFO("HapCol");
#endif
  INFO("Starting...");
  Counter MAX_COV = 0;
  Counter coverage = 0;
  Counter MAX_L = 0;
  Counter MAX_K = 0;
  Counter num_col = 0;
  vector<Counter> sum_successive_L;
  vector<vector<Counter> > scheme_backtrace;
  
  const constants_t constants;

  const options_t options= parse_arguments(argc, argv);
  INFO("Arguments:");
  INFO("Initialized? " << (options.options_initialized?"True":"False"));
  INFO("Input filename: '" << options.input_filename << '\'');
  INFO("Haplotype filename: '" << options.haplotype_filename << '\'');
  INFO("Discard weights? " << (options.unweighted?"True":"False"));
  INFO("Do not add X's? " << (options.no_xs?"True":"False"));
  INFO("Error rate: " << options.error_rate);
  INFO("Alpha: " << options.alpha);

  if (!options.options_initialized) {
    FATAL("Arguments not correctly initialized! Exiting..");
    exit(EXIT_FAILURE);
  }

  //.:: COLLECT STARTING PARAMETERS

  //Initializing the starting parameters: no competitive section
  
  //Pre-compute binomial values
  binom_coeff::initialize_binomial_coefficients(MAX_COVERAGE, MAX_COVERAGE);
  computeK(MAX_COVERAGE, options.alpha, options.error_rate);
  
  computeInputParams(num_col, MAX_COV, MAX_L, MAX_K, sum_successive_L,
                     coverage, options.input_filename, options.unweighted, scheme_backtrace);

  DEBUG(">> Initialized starting parameters");
  INFO("::== Starting parameters:  MAX_COV = " << MAX_COV << " // MAX_L = " << MAX_L << " // MAX_K = " << MAX_K);
  INFO("::== no of columns:     " << num_col);
  DEBUG("-->> sum_successive_L:  " << sum_successive_L);

  //.:: ALLOCATION MEMORY

  DEBUG(">> Starting allocation of memory");
  //Allocation of memory for input window
  vector<Column> input(2 * (MAX_L - 1) + 1,
                       Column(MAX_COV,
                              Entry(-1, Entry::BLANK, -1)));
  Pointer input_pointer = 0;
  TRACE("-->> input allocated");


  //Allocation of memory for backward indexer
  //The index in j of the shared elements between p and j
  vector<vector<Pointer> > back_indexer(2 * (MAX_L - 1) + 1,
                                      vector<Pointer>(MAX_COV,
                                                      -1));
  //Equal to indexer_pointer
  TRACE("-->> back indexer allocated");

  //Allocation of memory for forward indexer
  //The index in p of the shared elements between p and j
  vector<vector<Pointer> > forw_indexer(2 * (MAX_L - 1) + 1,
                                      vector<Pointer>(MAX_COV,
                                                      -1));
  const Pointer indexer_pointer = MAX_L - 1;
  TRACE("-->> forw indexer allocated");

  //Allocation of memory for vector of k_j
  vector<Counter> k_j(2 * (MAX_L - 1) + 1,
                      MAX_K);
  //its pointer is equal to input_pointer
  TRACE("-->> k_j allocated");

  //Allocation of memory for homozigous costs
  vector<Counter> homo_cost(2 * (MAX_L - 1) + 1,
                            MAX_COUNTER);
  //its pointer is equal to input_pointer
  TRACE("-->> homo_cost allocated");

  //Allocation of memory for homozigous weights
  vector<Cost> homo_weight(2 * (MAX_L - 1) + 1,
                           Cost::INFTY);
  //its pointer is equal to input_pointer
  TRACE("-->> homo_weight allocated");

  //Allocation of memory for prevision matrix
  //[Destinatary of prevision][Who make the prevision][Indexof(mask of who makes prevision on common fragments)]

  vector<vector<vector<Cost> > > prevision(MAX_L,
                                           vector<vector<Cost> > (MAX_L,
                                                                  vector<Cost>(0)));
  
  for(unsigned int j = 0; j < MAX_L; j++) {
    for(unsigned int q = 0; q < MAX_L; q++) {
      prevision[j][q].resize(sum_successive_L[q],
                             Cost::INFTY);
    }
  }
  Pointer prevision_pointer = 0;
  TRACE("-->> prevision allocated");

  //Allocation of memory for OPT vector
  vector<Cost> OPT(MAX_L + 1,
                   Cost::INFTY);        //+ 1 since I need OPT[j - L]
  Pointer OPT_pointer = 0;
  TRACE("-->> OPT allocated");

  //Allocation of memory for backtrace column
  //XXX: IMPROVEMENT: I can compute the precise sum_successive_L and MAX_L for each column of the input
  vector<vector<vector<Backtrace1> > > backtrace_table1(num_col);
  vector<vector<vector<bool> > > backtrace_table2_haplotypes(num_col);
  vector<vector<vector<bool> > > backtrace_table2_new_block(num_col);

  for(unsigned int j = 0; j < num_col; j++) {
    backtrace_table1[j].resize(scheme_backtrace[j].size());
    backtrace_table2_haplotypes[j].resize(scheme_backtrace[j].size());
    backtrace_table2_new_block[j].resize(scheme_backtrace[j].size());
    for(unsigned int q = 0; q < backtrace_table1[j].size(); q++) {
      backtrace_table1[j][q].resize(scheme_backtrace[j][q]);
      backtrace_table2_haplotypes[j][q].resize(scheme_backtrace[j][q]);
      backtrace_table2_new_block[j][q].resize(scheme_backtrace[j][q]);
    }
  }
  TRACE("-->> Backtrace table allocated");

  vector<bool> is_homozygous(num_col);
  vector<bool> homo_haplotypes(num_col);
  vector<Backtrace1> best_heterozygous1(num_col);
  vector<bool> best_heterozygous2_haplotypes(num_col);
  vector<bool> best_heterozygous2_new_block(num_col);

  
  DEBUG(">> Completed allocation of memory");

  //INITIALIZATION

  Combinations generator;
  ColumnReader column_reader(options.input_filename, coverage, options.unweighted);
  Column column;
  Counter step = 0;
  
  
  //Place the first and the next L columns in the positions of input data structure

  input_pointer = 0;

  Counter l = 0;

  while(column_reader.has_next() && l < MAX_L) {
    Pointer new_l_pointer = next(input_pointer, input.size(), l);

    if(l == 0) {
      column.clear();
    } else {
      column = column_reader.get_next();
    }

    insert_col_and_update(input, k_j, homo_cost, homo_weight, new_l_pointer, column, 
                          options.unweighted, homo_haplotypes, step + l);

    l++;
  }


  DEBUG(">> Initialization completed");


  //  .::: BASE CASE :::.
  INFO(".:: Basic Step: " << step);

  BitColumn colj;
  BitColumn corrected_colj;
  BitColumn mask;
  Cost current_cost;
  Cost current_best(Cost::INFTY);
  Counter cov_j;
  bool feasibility;
  bool has_successive;
  bool solution_existence(true);
  Counter temp_jump(-1);
  Counter temp_index(0);
  bool temp_haplotypes(false);
  bool temp_new_block(false);

  //Base case for OPT
  OPT[OPT_pointer] = 0;

  //k_j, homo_weight and homo_cost for base case
  k_j[input_pointer] = 0;
  homo_weight[input_pointer] = 0;
  homo_cost[input_pointer] = 0;

  //Initialize D[j, C'j] to infinite
  current_cost = 0;

  //Make a prevision for all the successive column
  has_successive = true;
  Counter p = 1;

  do {
    //All the condition that have to be satisfied for the next column
    const Pointer new_p_pointer = next(input_pointer, input.size(), p - 1);
    feasibility = (p - 1 == 0) || (homo_cost[new_p_pointer] <= k_j[new_p_pointer]);

    if (p >= MAX_L || forw_indexer[indexer_pointer + p][0] == -1 || !feasibility) {
      has_successive = false;
    } else {
      //The number of elements shared between p and j
      Pointer new_prevision_pointer = next(prevision_pointer, prevision.size(), p);
            
      prevision[new_prevision_pointer][p].resize(1);
      prevision[new_prevision_pointer][p][0] = current_cost;
      
      p++;
    }
  } while (has_successive);
  

  DEBUG("-->> Basic case completed  -- current_cost: " << current_cost);
  if (step % 500 == 0) {
    INFO(".:: Step: " << step << "  ==>  OPT: " << OPT[OPT_pointer]);
  } else {
    DEBUG(".:: Step: " << step << "  ==>  OPT: " << OPT[OPT_pointer]);
  }


  //DP

  //For all the columns
  while(!check_end(column_reader, input, next(input_pointer, input.size(), 1)) && solution_existence)
    {
      current_best = Cost::INFTY;
      solution_existence = false;
      temp_jump = -1;
      temp_index = 0;
      temp_haplotypes = false;
      temp_new_block = false;
      step++;

      // >>>>>>>>>>>>>>>>>>>>>> UPDATE DATA STRUCTURE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      //.:: Read Column
      if(column_reader.has_next()) {
        column = column_reader.get_next();
      } else {
        column.clear();
      }

      //.:: Update input

      //.:: Update common pointer
      input_pointer = next(input_pointer, input.size(), 1);

      Pointer new_input_pointer = next(input_pointer, input.size(), MAX_L - 1);

      insert_col_and_update(input, k_j, homo_cost, homo_weight, new_input_pointer, 
                            column, options.unweighted, homo_haplotypes, step + (MAX_L - 1));

      //.:: Update indexers

      //For all the q successive columns
      for(unsigned int q = 1; q < MAX_L; q++)
        {
          intersect(input[next(input_pointer, input.size(), q)],
                    input[input_pointer],
                    indexer_pointer + q,
                    forw_indexer, back_indexer);

          //If the just considered column q did not have any common element, we can fill as
          //  empty all the remaining columns and terminate.
          if(forw_indexer[indexer_pointer + q][0] == -1)
            {
              for(unsigned int p = q + 1; p < MAX_L; p++)
                {
                  forw_indexer[indexer_pointer + p][0] = -1;
                  back_indexer[indexer_pointer + p][0] = -1;
                }
              q = MAX_L;
            }
        }

      //For all the previous q columns
      for(unsigned int q = 1; q < MAX_L; q++)
        {
          intersect(input[prev(input_pointer, input.size(), q)],
                    input[input_pointer],
                    indexer_pointer - q,
                    forw_indexer, back_indexer);

          //If the just considered column q did not have any common element, we can fill as
          //  empty all the remaining columns and terminate.
          if(forw_indexer[indexer_pointer -  q][0] == -1)
            {
              for(unsigned int p = q + 1; p < MAX_L; p++)
                {
                  forw_indexer[indexer_pointer - p][0] = -1;
                  back_indexer[indexer_pointer - p][0] = -1;
                }
              q = MAX_L;
            }
        }


      //.:: Update prevision
      //XXX: Ricontrollare
      prevision_pointer = next(prevision_pointer, prevision.size(), 1);
      Pointer new_prevision_pointer = next(prevision_pointer, prevision.size(), MAX_L - 1);
      Pointer last_input_pointer = next(input_pointer, input.size(), MAX_L - 1);

      //Notice: prevision[new_prevision_pointer].size() == MAX_L
      for(unsigned int i = 1; i < prevision[new_prevision_pointer].size(); i++)
        {
          Pointer prec_input_pointer = prev(last_input_pointer, input.size(), i);
          Counter active_common = compute_active_common(input[prec_input_pointer], input[last_input_pointer]);
          
          /**
          TRACE("-----TEMP-----  active common: " << active_common << " k_j: " << k_j[prec_input_pointer]);
          TRACE("-----TEMP-----  total: " << cumulative_binomial_coefficient(active_common,
                                                                            k_j[prec_input_pointer])
               );
          TRACE("-----TEMP-----  but....: " << prevision[new_prevision_pointer][i].size());
          **/

          fill(prevision[new_prevision_pointer][i].begin(),
               (prevision[new_prevision_pointer][i].begin() +
                binom_coeff::cumulative_binomial_coefficient(active_common,
                                                             k_j[prec_input_pointer])),
               Cost::INFTY);
        }

      //.:: Update OPT

      OPT_pointer = next(OPT_pointer, OPT.size(), 1);
      OPT[OPT_pointer] = Cost::INFTY;

      
      DEBUG(">> Update data structure completed");

      //>>>>>>>>>>>>>>>>>>>>> END UPDATE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      //>>>>>>>>>>>>>>>>>>>>> ITERATIVE STEP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      //Binary repesentation of the column and compute of coverage
      represent_column(input[input_pointer], colj, cov_j);

      DEBUG("...| Column: " <<  column_to_string(colj, cov_j) << " -- current coverage: " << cov_j << " and current k: " << k_j[input_pointer]);
      
      //Initializing OPT[j] = infinite
      //XXX: Is it redundant??
      //OPT[OPT_pointer] = Cost::INFTY;

      //We have already computed k, homo_weights and homo_cost for current column
      
      //First option for the value of OPT[j]
      if(homo_cost[input_pointer] <= k_j[input_pointer])
        {
          //XXX: Can I remove this check?
          Cost temp = homo_weight[input_pointer] + OPT[prev(OPT_pointer, OPT.size(), 1)];
          if(temp < OPT[OPT_pointer]) {
            OPT[OPT_pointer] =  temp;
            solution_existence = true;
            is_homozygous[step] = true;
            DEBUG(".:: Column: " << step << " can be homozygous with a cost: " << OPT[OPT_pointer]);
          }
        }


      //Enumerate all the combinations

      generator.initialize_cumulative(cov_j, k_j[input_pointer]);
      while(generator.has_next())
        {
          generator.next();
          generator.get_combination(mask);
          TRACE("|--------");
          TRACE("|== Mask: " << column_to_string(mask, cov_j));

          //Initialize D[j, C'j] to infinite
          current_cost = Cost::INFTY;

          //Compute C'j
          corrected_colj = colj ^ mask;
          TRACE("-->> corrected column: " << column_to_string(corrected_colj, cov_j));

          //The column cannot be transformed into an homozygous column
          if(corrected_colj.any() && (corrected_colj.count() != cov_j) )
            {
              //Compute the weight of the mask
              Cost weight_mask = 0;

              if (options.unweighted) {
                weight_mask = Cost((Cost::cost_t)mask.count());
              } else {
                compute_weight_mask(mask, input[input_pointer], weight_mask);
              }

              //Compute current_cost that corresponds to D[j, Bj]

              Counter q = 1;
              Pointer new_homo_pointer = prev(input_pointer, input.size(), q - 1);
              bool has_previous = true;
              Cost cumulative_homo = 0;

              do {
                //All the condition that have to be satisfied for the next column
                feasibility = (q - 1 == 0) || (homo_cost[new_homo_pointer] <= k_j[new_homo_pointer]);

                if (q >= MAX_L || forw_indexer[indexer_pointer - q][0] == -1 || !feasibility) {
                  has_previous = false;
                } else {
                  Counter active_qj(0);
                  BitColumn mask_qj;
                  Cost temp(0);

                  Pointer new_q_pointer = prev(input_pointer, input.size(), q);

                  //First Mask
                  extract_common_mask(input[new_q_pointer], indexer_pointer - q, input[input_pointer],
                                      mask, back_indexer, forw_indexer, mask_qj, active_qj);

                  if(mask_qj.count() <= k_j[new_q_pointer])
                    {
                      Counter index = generator.cumulative_indexof(mask_qj, active_qj);
                      temp = prevision[prevision_pointer][q][index] + weight_mask + cumulative_homo;
                      if(temp < current_cost) {
                        current_cost = temp;
                        solution_existence = true;

                        temp_jump = q;
                        temp_index = index;
                        temp_haplotypes = backtrace_table2_haplotypes[step - q][q][index];
                        temp_new_block = false;
                      }
                      TRACE("-->> Temporary current cost: " << current_cost); 
                      TRACE("---->> the previous equal heterozigous is " << (step - q) 
                            << "  -- its mask: " << column_to_string(mask_qj, active_qj));
                    }

                  //Complement
                  complement_mask(mask_qj, active_qj, constants);

                  if(mask_qj.count() <= k_j[new_q_pointer])
                    {
                      Counter index = generator.cumulative_indexof(mask_qj, active_qj);
                      temp = prevision[prevision_pointer][q][index] + weight_mask + cumulative_homo;
                      if(temp < current_cost) {
                        current_cost = temp;
                        solution_existence = true;
                        
                        temp_jump = q;
                        temp_index = index;
                        temp_haplotypes = !backtrace_table2_haplotypes[step - q][q][index];
                        temp_new_block = false;
                      }
                      TRACE("-->> Temporary current cost: " << current_cost); 
                      TRACE("---->> the previous equal heterozigous is " << (step - q) 
                            << "  -- its mask: " << column_to_string(mask_qj, active_qj));
                    }

                  q++;

                  new_homo_pointer = prev(input_pointer, input.size(), q - 1);
                  
                  cumulative_homo += homo_weight[new_homo_pointer];
                  TRACE("----> Cumulative homo: " << cumulative_homo << "  with q:  " << q);
                }
              } while(has_previous);

              TRACE("-->> Best current cost (D[j, C'j]): "<< current_cost);

              //Third case of the recursion for D[j, C'j]
              //XXX: Check carefully!
              if(q <= MAX_L && feasibility) {
                Cost temp = OPT[prev(OPT_pointer, OPT.size(), q)] + weight_mask + cumulative_homo;
                if(temp < current_cost) {
                  current_cost = temp;
                  solution_existence = true;

                  temp_jump = q;
                  temp_index = 0;
                  temp_haplotypes = false;
                  temp_new_block = true;
                  TRACE("<<>> Third case of recursion - First heterozigous of new block");
                  TRACE("..OPT[previous] = " << OPT[prev(OPT_pointer, OPT.size(), q)] << " - weight:  " << weight_mask << " - cumulative_homo: " << cumulative_homo);
                  //TRACE(".:: Column: " << step << " can be heterozigous with a cost: " << current_cost);
                  //TRACE("====> Best correction:  " << column_to_string(mask, cov_j));
                }
              }

              //Make a prevision for all the seccessive column
              has_successive = true;
              Counter p = 1;

              //All the condition that have to be satisfied for the next column

              do {
                Pointer new_homo_pointer = next(input_pointer, input.size(), p - 1);
                feasibility = (p - 1 == 0) || (homo_cost[new_homo_pointer] <= k_j[new_homo_pointer]);

                if (p >= MAX_L || forw_indexer[indexer_pointer + p][0] == -1 || !feasibility) {
                  has_successive = false;
                } else {
                  //The number of elements shared between p and j
                  Counter active_pj = 0;
                  cut(mask, back_indexer[indexer_pointer + p], active_pj);
                  TRACE("-->> Successive column: " << (step + p) << " -- Common elements:  " << active_pj << " -- Cut mask: " << column_to_string(mask, active_pj));

                  Counter index = generator.cumulative_indexof(mask, active_pj);
                  Pointer new_prevision_pointer = next(prevision_pointer, prevision.size(), p);
                  Cost& temp = prevision[new_prevision_pointer][p][index];
                  if(current_cost < temp) {
                    temp = current_cost;

                    backtrace_table1[step][p][index].jump = temp_jump;
                    backtrace_table1[step][p][index].index = temp_index;
                    backtrace_table2_haplotypes[step][p][index] = temp_haplotypes;
                    backtrace_table2_new_block[step][p][index] = temp_new_block;
                  }
                  TRACE("USCITO:  ");
                  p++;
                }
              } while(has_successive);

              
              if(current_cost < current_best) {
                current_best = current_cost;

                best_heterozygous1[step].jump = temp_jump;
                best_heterozygous1[step].index = temp_index;
                best_heterozygous2_haplotypes[step] = temp_haplotypes;
                best_heterozygous2_new_block[step] = temp_new_block;
              }
              

              //Update value of OPT for the current column
              if(current_cost < OPT[OPT_pointer]) {
                OPT[OPT_pointer] = current_cost;
                is_homozygous[step] = false;
                DEBUG(".:: Column: " << step << " can be heterozigous with a cost: " << OPT[OPT_pointer]);
                DEBUG("====> Best correction:  " << column_to_string(mask, cov_j));
              }
              TRACE("-->> OPT: " << OPT[OPT_pointer]);
            }
        }

      if (step % 500 == 0) {
        INFO(".:: Step: " << step << "  ==>  OPT: " << OPT[OPT_pointer]);
      } else {
        DEBUG(".:: Step: " << step << "  ==>  OPT: " << OPT[OPT_pointer]);
      }
      //End of DP cycle for all the columns
    }
  
  if(solution_existence) {
    INFO("*** SUCCESS ***");
    INFO("===> Optimal cost:  " << OPT[OPT_pointer]);
    vector<bool> haplotype1(num_col - 1);
    vector<bool> haplotype2(num_col - 1);
    reconstruct_haplotypes(backtrace_table1, backtrace_table2_haplotypes, backtrace_table2_new_block,
                           is_homozygous, homo_haplotypes,
                           best_heterozygous1, best_heterozygous2_haplotypes, best_heterozygous2_new_block,
                           haplotype1, haplotype2);

    //INFO("<=> Haplotype 1:  " << haplotype_1);
    //INFO("<=> Haplotype 2:  " << haplotype_2);

    vector<char> output1(haplotype1.size());
    vector<char> output2(haplotype2.size());
    

    if(!options.no_xs) {
      add_xs(haplotype1, haplotype2, output1, output2, options.input_filename, coverage, options.unweighted);
    }

    DEBUG("<<>> Writing haplotypes...");
		ofstream ofs;
		try { 
			ofs.open(options.haplotype_filename.c_str(), ios::out);
      if(!options.no_xs) {
        ofs << output1 << endl;
        ofs << output2 << endl;
      } else {
        ofs << haplotype1 << endl;
        ofs << haplotype2 << endl;
      }
    } catch(exception & e) { 
			ERROR("::::::: Error writing haplotype to \"" << options.haplotype_filename << "\": " << e.what());
			if(!options.no_xs) {
        ofs << output1 << endl;
        ofs << output2 << endl;
      } else {
        ofs << haplotype1 << endl;
        ofs << haplotype2 << endl;
      }
      return EXIT_FAILURE;
		}
  } else {
    INFO("*** NO SOLUTION ***");
    INFO("<<>> No feasible solution exist with these parameters -- alpha = " << options.alpha << " and error rate = " << options.error_rate); 
    INFO("<<>> The last not feasible column is:  " << step << "  with coverage = " << cov_j << " and k = " << k_j[input_pointer]);
  }
}




//Change such that we do not take all the input in memory!!
void computeInputParams(Counter &num_cols, Counter &MAX_COV, Counter &MAX_L,
                        Counter &MAX_K, vector<Counter> &sum_successive_L,
                        Counter &coverage, const string &inputfilename, const bool &unweighted,
                        vector<vector<Counter> > &scheme_backtrace)
{

  //XXX: Not MAX_COVERAGE but coverage, so we can choose the coverage that we want to obtain
  coverage = MAX_COVERAGE;

  ColumnReader column_reader(inputfilename, coverage, unweighted);

  num_cols = column_reader.num_cols() + 1; //We add a starting dummy empty column

  vector<Column> input(num_cols);
  vector<unsigned int> homo_cost(num_cols);
  vector<Column>::iterator input_iterator(input.begin());
  vector<Counter> rows(num_cols * MAX_COVERAGE, 0);

  while(input_iterator != input.end() && column_reader.has_next())
    {
      Counter count_major = 0;
      Counter count_minor = 0;
      
      Column &current_column = *input_iterator;
      Column read_column;

      if(input_iterator == input.begin()) {
        read_column.clear();
      } else {
        read_column = column_reader.get_next();
      }

      current_column.resize(read_column.size(), Entry(-1, Entry::BLANK, -1));
      
      for(unsigned int i = 0; i < read_column.size(); ++i) {
        current_column[i].set_read_id(read_column[i].get_read_id());
        current_column[i].set_allele_type(read_column[i].get_allele_type());
        current_column[i].set_phred_score(read_column[i].get_phred_score());

        if(current_column[i].get_allele_type() == Entry::MAJOR_ALLELE) {
          ++count_major;
        } else {
          ++count_minor;
        }

        ++rows[current_column[i].get_read_id()];
      }

      //sufficient condition to check the feasibility for the homozygous transformation
      homo_cost[input_iterator - input.begin()] = std::min(count_major, count_minor);

      MAX_COV = std::max(static_cast<Counter>((*input_iterator).size()), MAX_COV);

      ++input_iterator;
    }

  MAX_L = *max_element(rows.begin(), rows.end());
  MAX_K = computeK(MAX_COV);
  Counter MAX_CONS_HOMO = 0;    //The maximum number of consecutive homozigous columns
  
  sum_successive_L.resize(MAX_L, 0);
  scheme_backtrace.clear();
  scheme_backtrace.resize(num_cols);
  for(unsigned int i = 0; i < input.size(); i++)
    {
      unsigned int y = 1;
      unsigned int k_temp = computeK(input[i].size());
      unsigned int current_cons_homo = 0;   //The maximum number of consecutive homozigous columns assuming i the first
      bool flag = true;
      scheme_backtrace[i].push_back(0);

      while(y < MAX_L && (i + y) < input.size())
        {
          Counter active_common = compute_active_common(input[i], input[i + y]);

          Counter result = binom_coeff::cumulative_binomial_coefficient(active_common,
                                                                        k_temp);
          sum_successive_L[y] = max(sum_successive_L[y], result);
          
          if(flag) {
            //XXX: Can I add && active_common != 0?
            if( (homo_cost[i + y] <= computeK(input[i + y].size())) && active_common != 0) {
              ++current_cons_homo;
              scheme_backtrace[i].push_back(result);
            } else {
              flag = false;
              scheme_backtrace[i].push_back(result);
            }
          }

          y++;
        }
      
      MAX_CONS_HOMO = std::max(MAX_CONS_HOMO, current_cons_homo);
    }

  //+1 is necessary to count the first heterozygous column before the longest sequence of homozugouses
  //and another +1 to count the heterozygous column after that

  //min it is necessary since MAX_CONS_HOMO + 2 can be larger than the real MAX_L
  MAX_L = std::min(MAX_L, MAX_CONS_HOMO + 2);
}


//XXX: Can I leave parameter q and pass as parameter just one column of forw_indexer and back_indexer??????????
void intersect(const Column &colQ, const Column &colJ, const Pointer &q,
               vector<vector<Pointer> > &forw_indexer, vector<vector<Pointer> > &back_indexer)
{
  size_t i = 0;
  size_t j = 0;
  size_t count = 0;

  while ((i < colQ.size()) && 
         (j < colJ.size()) &&
         (colJ[j].get_read_id() != -1) && 
         (colQ[i].get_read_id() != -1)) {
    if (colQ[i].get_read_id() == colJ[j].get_read_id()) {
      forw_indexer[q][count] = i;
      back_indexer[q][count] = j;
      ++i;
      ++j;
      ++count;
    } else if (colQ[i].get_read_id() < colJ[j].get_read_id()) {
      ++i;
    } else {
      ++j;
    }
  }
  
  if(count < forw_indexer[q].size()) {
    forw_indexer[q][count] = -1;
    back_indexer[q][count] = -1;
  }
}


void represent_column(const Column &column, BitColumn &result, Counter &cov)
{
  result.reset();
  cov = 0;
  while(cov < column.size() && column[cov].get_read_id() != -1) {
    result.set(cov, (column[cov].get_allele_type() == Entry::MINOR_ALLELE));
    ++cov;
  }
}


void cut(BitColumn &to_be_cut, const vector<Pointer> &indexer, Counter &active_pj)
{
  active_pj = 0;
  const BitColumn copy(to_be_cut);

  to_be_cut.reset();

  const size_t isize= indexer.size();
  const Pointer* iactive= &indexer[0];
  while(active_pj < isize && *iactive != -1) {
    to_be_cut.set(active_pj, copy[*iactive]);
    ++active_pj;
    ++iactive;
  }
}



void extract_common_mask(const Column &column_q, const Pointer &q_pointer,
                         const Column &column_j, const BitColumn &mask_colj,
                         const vector<vector<Pointer> > &back_indexer,
                         const vector<vector<Pointer> > &forw_indexer,
                         BitColumn &mask_qj, Counter &active_qj)
{
  mask_qj.reset();
  active_qj = 0;

  const Pointer* forw_indexer_q= &forw_indexer[q_pointer][0];
  const Pointer* back_indexer_q= &back_indexer[q_pointer][0];
  const size_t back_indexer_size= back_indexer[q_pointer].size();
  while (active_qj < back_indexer_size && *back_indexer_q != -1) {
    if ((column_q[*forw_indexer_q].get_allele_type() !=
         column_j[*back_indexer_q].get_allele_type())
        != mask_colj.test(*back_indexer_q)) {
      mask_qj.set(active_qj, 1);
    }
    ++active_qj;
    ++forw_indexer_q;
    ++back_indexer_q;
  }
}


int compute_active_common(const Column &colJ, const Column &colQ)
{
  int active_common = 0;
  unsigned int i = 0;
  unsigned int j = 0;
  
  while ((i < colQ.size()) && 
         (j < colJ.size()) &&
         (colJ[j].get_read_id() != -1) && 
         (colQ[i].get_read_id() != -1)) {
    if (colQ[i].get_read_id() == colJ[j].get_read_id()) {
      ++i;
      ++j;
      ++active_common;
    } else if (colQ[i].get_read_id() < colJ[j].get_read_id()) {
      ++i;
    } else {
      ++j;
    }
  }
  return active_common;
}


void insert_col_and_update(vector<Column> &input, vector<Counter> &k_j, vector <Counter> &homo_cost,
                           vector<Cost> &homo_weight, const Pointer &pointer, const Column &column, const bool &unweighted,
                           vector<bool> &kind_homozygous, const Counter &step) {
  
  Counter count_major = 0;
  Cost weight_major = 0;

  Counter count_minor = 0;
  Cost weight_minor = 0;

  Counter column_size = column.size();

  unsigned int i = 0;
  for(i = 0; i < column.size(); i++)
    {
      const readid_t column_read_id = column[i].get_read_id();
      const Entry::allele_t column_allele_type = column[i].get_allele_type();
      const unsigned int column_phred_score = (unweighted)? 1 : column[i].get_phred_score();

      input[pointer][i].set_read_id(column_read_id);
      input[pointer][i].set_allele_type(column_allele_type);
      input[pointer][i].set_phred_score(column_phred_score);

      if(column_allele_type == Entry::MINOR_ALLELE) {
        count_minor++;
        weight_minor += column_phred_score;
      } else if (column_allele_type == Entry::MAJOR_ALLELE) {
        count_major++;
        weight_major += column_phred_score;
      } else {
        throw runtime_error("the input data contains an allele that is not equal to 0 or 1");
      }
    }

  if(i < input[pointer].size() && input[pointer][i].get_read_id() != -1)
    {
      input[pointer][i].set_read_id(-1);
      input[pointer][i].set_allele_type(Entry::BLANK);
      input[pointer][i].set_phred_score(-1);
    }

  //.:: Update k_j for current column

  k_j[pointer] = computeK(column_size);


  //.:: Update homozygous cost
      
  homo_cost[pointer] = MAX_COUNTER;
  homo_weight[pointer] = Cost::INFTY;

  if(count_minor <= k_j[pointer] && weight_minor < homo_weight[pointer]) {
    homo_cost[pointer] = count_minor;
    homo_weight[pointer] = weight_minor;
    if(step < kind_homozygous.size()) {
      kind_homozygous[step] = true;
    }
  }

  if(count_major <= k_j[pointer] && weight_major < homo_weight[pointer]) {
    homo_cost[pointer] = count_major;
    homo_weight[pointer] = weight_major;
    if(step < kind_homozygous.size()) {
      kind_homozygous[step] = false;
    }
  }
}



void compute_weight_mask(const BitColumn &mask, const Column &column, Cost &weight_mask) {
  BitColumn comb(mask);
  weight_mask = 0;  
  int pos = 0;
  int temp = 0;
    
  while(comb.any())
    {
      temp = ffsl(comb.to_ulong());
      pos += temp;
      weight_mask += column[pos - 1].get_phred_score();
      comb>>=(temp);
    }
}


void reconstruct_haplotypes(const vector<vector<vector<Backtrace1> > > &backtrace_table1,
                            const vector<vector<vector<bool> > > &backtrace_table2_haplotypes,
                            const vector<vector<vector<bool> > > &backtrace_table2_new_block,
                            const vector<bool> &is_homozygous,
                            const vector<bool> &homo_haplotypes,
                            const vector<Backtrace1> &best_heterozygous1,
                            const vector<bool> &best_heterozygous2_haplotypes,
                            const vector<bool> &best_heterozygous2_new_block,
                            vector<bool> &haplotype1, vector<bool> &haplotype2) {
  Counter col = backtrace_table1.size() - 1;

  haplotype1.resize(col);
  haplotype2.resize(col);

  while(col > 0) {
    while(is_homozygous[col]) {
      if(homo_haplotypes[col]) {
        haplotype1[col - 1] = false;
        haplotype2[col - 1] = false;
      } else {
        haplotype1[col - 1] = true;
        haplotype2[col - 1] = true;
      }

      --col;
    }
  
    Backtrace1 back1 = best_heterozygous1[col];
    bool back2_haplotypes = best_heterozygous2_haplotypes[col];
    bool back2_new_block = best_heterozygous2_new_block[col];
    bool flag = col > 0;

    while (flag) {
      if(back2_haplotypes) {
        haplotype1[col - 1] = false;
        haplotype2[col - 1] = true;
      } else {
        haplotype1[col - 1] = true;
        haplotype2[col - 1] = false;
      }

      for(int i = 0; i < (back1.jump - 1); i++) {
        --col;
        if(homo_haplotypes[col]) {
          haplotype1[col - 1] = false;
          haplotype2[col - 1] = false;
        } else {
          haplotype1[col - 1] = true;
          haplotype2[col - 1] = true;
        }
      }
      
      --col;
      
      if(back2_new_block || col == 0) {
        flag = false;
      } else {
        flag = true;
        back2_haplotypes = backtrace_table2_haplotypes[col][back1.jump][back1.index];
        back2_new_block = backtrace_table2_new_block[col][back1.jump][back1.index];
        back1 = backtrace_table1[col][back1.jump][back1.index];
      }
      
    }
  }
}



Counter computeK(const Counter &cov, const double &alpha, const double &error_rate)
{
  static bool computed = false;
  static vector<Counter> ks(cov + 1, 0);

  if (!computed) {
    for(Counter i = 1; i < ks.size(); ++i) {
      Counter k = 0;

      double cumulative =  pow(1.0 - error_rate, i);

      while(!(1.0 - cumulative <= alpha) && (k < i)) {
        ++k;
        cumulative += (double)binom_coeff::binomial_coefficient(i, k) * pow(error_rate, k) * pow(1.0 - error_rate, i - k);
      }
      
      ks[i] = k;
    }
    computed = true;
  }
    
  return ks[cov];
}



void add_xs(const vector<bool> &haplo1, const vector<bool> &haplo2, 
            vector<char> &haplo1_out, vector<char> &haplo2_out,
            const string &input_filename, const unsigned int coverage, const bool unweighted) {
  
  ColumnReader column_reader(input_filename, coverage, unweighted);

  vector<vector<bool> > reads_matrix;
  vector<vector<unsigned int> > weights;
  vector<int> starting_positions;

  vector<vector<bool> > mapping_haplo1(column_reader.num_cols());
  vector<vector<bool> > mapping_haplo2(column_reader.num_cols());

  int current_column = -1;

  while(column_reader.has_next()) {
    vector<Entry> column = column_reader.get_next();
    
    ++current_column;

    for(unsigned int i = 0; i < column.size(); ++i) {
      while(reads_matrix.size() <= (unsigned int)column[i].get_read_id()) {
        reads_matrix.push_back(vector<bool>(0));
        weights.push_back(vector<unsigned int>(0));
        starting_positions.push_back(-1);
      }

      if(reads_matrix[column[i].get_read_id()].size() == 0) {
        starting_positions[column[i].get_read_id()] = current_column;
      }

      reads_matrix[column[i].get_read_id()].push_back(column[i].get_allele_type() == Entry::MINOR_ALLELE);
      
      if(unweighted) {
        weights[column[i].get_read_id()].push_back(1);
      } else {
        weights[column[i].get_read_id()].push_back(column[i].get_phred_score());
      }
    }

  }
  
  unsigned int total_errors = 0;


  for(unsigned int read_id = 0; read_id < reads_matrix.size(); ++read_id) {
    if(map_fragment(reads_matrix[read_id], weights[read_id], starting_positions[read_id], haplo1, haplo2, total_errors) == 1) {
      
      //Map read in haplo 1
      for(unsigned col = 0; col < reads_matrix[read_id].size(); ++col) {
        mapping_haplo1[col + starting_positions[read_id]].push_back(reads_matrix[read_id][col]);
      }
      
    } else {
      
      //Map read in haplo 2
      for(unsigned col = 0; col < reads_matrix[read_id].size(); ++col) {
        mapping_haplo2[col + starting_positions[read_id]].push_back(reads_matrix[read_id][col]);
      }
      
    }
  }

  make_haplo(haplo1, mapping_haplo1, haplo1_out);
  make_haplo(haplo2, mapping_haplo2, haplo2_out);

  INFO("TOTAL MISMATCHES DURING MAPPING:   " << total_errors);
}



int map_fragment(const vector<bool> &read, const vector<unsigned int> &weights, const unsigned int offset,
                 const vector<bool> &haplo1_in, const vector<bool> &haplo2_in, unsigned int &total_errors) {
  unsigned int distance1 = 0;
  unsigned int distance2 = 0;

  for(unsigned int col = 0; col < read.size(); ++col) {
    if(read[col] != haplo1_in[col + offset]) {
      distance1 += weights[col];
    } 
    if(read[col] != haplo2_in[col + offset]) {
      distance2 += weights[col];
    }    
  }

  if(distance1 <= distance2) {
    total_errors += distance1;
    return 1;
  } else {
    total_errors += distance2;
    return 2;
  }
}


void make_haplo(const vector<bool> &haplo, const vector<vector<bool> > &mapping_haplo, vector<char> &haplo_out) {
  unsigned int count_X = 0;
  vector<int> counter(2);

  for(unsigned int col = 0; col < mapping_haplo.size(); ++col) {
    count_alleles(mapping_haplo[col], counter);
    
    if(counter[0] == counter[1]) {
      haplo_out[col] = 'X';
      ++count_X;
    } else {
      if(haplo[col]) {
        haplo_out[col] = '1';
      } else {
        haplo_out[col] = '0';
      }
    }
  }

  INFO("INTRODUCED X's IN ONE HAPLOTYPE:   " << count_X);
}


void count_alleles(const vector<bool> &col, vector<int> &counter) {
  counter[0] = 0;
  counter[1] = 0;

  for(unsigned int i = 0; i < col.size(); ++i) {
    if(col[i] == false) {
      ++counter[0];
    } else if(col[i] == true) {
      ++counter[1];
    }
  }
}


