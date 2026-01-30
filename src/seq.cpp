#include "codon.h"
#include "seq.h"
#include <plog/Log.h>


codon::Seq::Seq(const std::string &input){

    int remainder_size      = input.length() % 3;
    bool all_codons_full    = (remainder_size == 0);

    this->seq.reserve((all_codons_full) ? input.length() / 3 : input.length() / 3 + 1);

PLOGD << "Generating Seq '" << input;
PLOGD << "Remainder_size = " << remainder_size << "\t//\tAll_codons_full = " << all_codons_full;


    for (int i = 0; i < input.length() / 3; i++) {
        seq.emplace_back(codon::Codon(input.substr(i*3, 3)));
        PLOGD << "Placing codon '" << input.substr(i*3, 3) << "' at pos. " << i;
    }

    if (!all_codons_full) {
        seq.emplace_back(codon::Codon(input.substr(input.length()-remainder_size)));
        PLOGD << "Placing incomplete codon '" << input.substr(input.length()-remainder_size) 
            << "' at pos. " << (input.length() / 3 + 1);
    }

}

codon::Seq::~Seq() {
    PLOGD << "Sequence at memory location '" << &this->seq << "' going out of scope";

}


std::string codon::Seq::get_seq_str() const {
    std::string annealed_str;
    annealed_str.reserve(this->seq.size()*3);

    for (codon::Codon curr_codon : this->seq) {
        annealed_str.append(curr_codon.get_bases_str());
        PLOGD << "Appending '" << curr_codon.get_bases_str() << "' to annealed_str.";
    }
    annealed_str.shrink_to_fit();
    PLOGD << "Generated final string: '" << annealed_str << "'.";
    return annealed_str;
};


void insert_base (codon::base  base,  std::size_t insert_loc, int shift_loc) {
//incomplete
}

void insert_codon (codon::Codon codon, std::size_t insert_loc, int shift_loc) {
//incomplete
}

std::size_t codon::Seq::get_seq_len() const {
    return this->seq.size();
}

void codon::Seq::insert_seq(codon::Seq other, std::size_t insert_loc, int shift_loc=0) {
    //check if capacity is compromised
    if ((this->seq.size() + other.seq.size()) < this->seq.capacity()) {
        this->seq.reserve((this->seq.size() + other.seq.size())*1.2);
        PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
    };
    int bases_at_insert = this->seq[insert_loc].get_bases_len();
    switch (bases_at_insert) {
        case 0:
            //replace the existing codon with the insert
        case 1: 
            if (shift_loc == 0) //extract base and insert at end of insert
            if (shift_loc == 1) //extract base and squeeze into start of insert
        case 2:
            if (shift_loc == 0) //extract bases 1 and 2 and insert at end of insert
            if (shift_loc == 1) //extract base 1 and insert at squeeze into start of insert
                                //extract base 2 and insert at the end of the insert
            if (shift_loc == 2) //extract base 1 and 2 and squeeze both into start of insert
        case 3:
            if (shift_loc == 0) //extract bases 1,2 and 3 and insert at end of insert
            if (shift_loc == 1) //extract base 1 and insert at squeeze into start of insert
                                //extract base 2 and 3 and insert at the end of the insert
            if (shift_loc == 2) //extract base 1 and 2 and squeeze both into start of insert
                                //extract base 3 and insert at the end of the insert
            if (shift_loc == 3) //extract base 1,2 and 3 and squeeze both into start of insert
    }
}


