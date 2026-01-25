#include "codon.h"
#include "seq.h"
#include <plog/Log.h>


codon::Seq::Seq(const std::string &input, int shift_start, int shift_stop){
    this->shift_startcodon  = shift_start;
    this->shift_stopcodon   = shift_stop;
    int remainder_size      = input.length() % 3;
    bool all_codons_full    = (remainder_size == 0);

    this->seq.reserve((all_codons_full) ? input.length() / 3 : input.length() / 3 + 1);
    auto it_start = input.begin();
    auto it_stop  = (input.length() > 3) ? it_start+3 : input.end();
    auto it_seq   = this->seq.begin();

    for (int i = 0; i < input.length() / 3; i++) {
        seq.emplace(it_seq++, codon::Codon(input.substr(i*3, 3)));
    }
    if (!all_codons_full) {
        seq.emplace(it_seq, codon::Codon(input.substr(input.length()-remainder_size)));
    }

}

codon::Seq::~Seq() {
    PLOGD << "Sequence at memory location '" << &this->seq << "' going out of scope";

}


std::string codon::Seq::get_encoding_str() const {
    std::string annealed_str;
    annealed_str.reserve(this->seq.size()*3);

    for (codon::Codon curr_codon : this->seq) {
        annealed_str.append(curr_codon.get_bases_str());
    }
    annealed_str.shrink_to_fit();
    return annealed_str;
};
