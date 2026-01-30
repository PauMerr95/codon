#pragma once
#include "codon.h"
#include "seq.h"
#include <string>
#include <vector>

namespace test {

    int codon_test();
    void check_creation_str(std::vector<std::string> arr_bases);
    void check_creation_str(std::vector<std::string> arr_bases, 
                            std::vector<codon::Codon> &generated);
    void check_creation_base(codon::base arr_bases[], int len);
    void check_operations(std::vector<codon::Codon> arr_codons);



    int seq_test();
    std::vector<codon::Seq> seq_build(const std::vector<std::string> &arr_sequences);

}
