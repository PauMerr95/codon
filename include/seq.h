#pragma once
#include <string>
#include "codon.h"
#include <vector>

namespace codon {

    class Seq {
        std::vector<codon::Codon> seq;

        public:
        Seq(const std::string &input);
        ~Seq();

        void insert_base  (codon::base base, std::size_t insert_loc, int shift_loc);
        void insert_codon (codon::Codon codon, std::size_t insert_loc, int shift_loc);
        void insert_seq   (codon::Seq other, std::size_t insert_loc, int shift_loc);
        
        void left_shift ();
        void right_shift ();

        std::string get_seq_str() const;
        codon::Seq  get_seq_base() const;
        std::vector<std::bitset<8>> get_seq_bin() const;




    };
}
