#include <string>
#include <algorithm>

#include "../src/balanced_combinations.h"


static inline
string column_to_string(const BitColumn &mask, const unsigned int &len) {
  string str = mask.to_string();
  reverse(str.begin(), str.end());
  return str.substr(0, len);
}


int main(int argc, char** argv)
{
    binom_coeff::initialize_binomial_coefficients(MAX_COVERAGE, MAX_COVERAGE);  
    BalancedCombinations generator;

    if(argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <COVERAGE> <COLUMN> <NUM_CORRECTIONS> <THRESHOLD>" << std::endl;
        std::cerr << "INFO:" << std::endl;
        std::cerr << "\t<COVERAGE>: an integer corresponding to the coverage" << std::endl;
        std::cerr << "\t<COLUMN>: an integer whose binary encoding corresponds to the starting column" << std::endl;
        std::cerr << "\t<NUM_CORRECTIONS>: an integer number of corrections to apply" << std::endl;
        std::cerr << "\t<THRESHOLD>: a double corresponding to the balancing threshold" << std::endl;

        return 1;
    }
    
    const unsigned int coverage(atoi(argv[1]));
    const unsigned int int_column(atoi(argv[2]));
    const bitset<32> column(atoi(argv[2]));
    const unsigned int corrections(atoi(argv[3]));
    const double threshold(atof(argv[4]));

    if(coverage < corrections)
    {
        std::cerr << "ERROR: the number of corrections must be at most the coverage!" << std::endl;
        return 0;
    }

    unsigned int check = int_column;
    unsigned int digits = 0;

    for (digits = 0; check > 0; check >>= 1)
        digits++;

    if(coverage < digits)
    {
        std::cerr << "ERROR: the coverage is not enough to represent a column as binary!" << std::endl;
        return 0;
    }

    if (threshold < 0.0 || threshold > 1.0)
    {
        std::cerr << "ERROR: the threshold must be within 0 and 1!" << std::endl;
        return 0;
    }

    std::cout << "Starting column:  " << column_to_string(column, coverage) << std::endl;
    
    generator.initialize(coverage, corrections, column, threshold);
    std::cout << "Results (in format CORRECTED_COLUMN::CORRECTIONS):  " << std::endl;
    while(generator.has_next())
    {
        bitset<32> result;
        generator.next();
        generator.get_combination(result);
        std::cout << "----------------  " << column_to_string(result^column, coverage) << " :: " << column_to_string(result, coverage)  << std::endl;
    }
}
