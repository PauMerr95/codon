#pragma once
#include <string>
#include "codon.h"
#include <vector>

namespace codon {

    class Seq {
        int shift_startcodon;
        int shift_stopcodon;
        std::vector<codon::Codon> seq;

        public:
        Seq(const std::string &input, int shift_start=0, int shift_stop=0);
        ~Seq();

        void insert_seq(codon::Seq other, int* insert_loc);
        void set_startcodon(int* startcodon_loc, int shift);
        void set_stopcodon(int* startcodon_loc, int shift);
        void realign();

        std::string get_encoding_str() const;
        codon::Seq  get_encoding_seq() const;
        std::vector<std::bitset<8>> get_encoding_bin() const;




    };
}
