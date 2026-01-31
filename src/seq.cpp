#include "codon.h"
#include "seq.h"
#include <plog/Log.h>


codon::Seq::Seq(const std::string &input){

    int remainder_size      = input.length() % 3;
    bool all_codons_full    = (remainder_size == 0);
    
    //ignore potential warning regarding integer division in combination with float -> int division part
    //calculates the amount of necessary codons and floating point multiplication is just for 20% safety buffer
    //loss of precision is negligable.
    this->seq.reserve(
        static_cast<std::size_t>(
            (all_codons_full) ? (input.length() / 3)*1.2 : (input.length() / 3 + 1)*1.2
        )
    );

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

void codon::Seq::left_shift() {
    //removes left most base in final codon and the starts squeeze chain to the front propogating
    //bases up until index 0 (only possible and meaningful) when first codon is incomplete).
    if (this->seq.at(0).get_bases_len() < 3) {
        std::size_t r_idx = this->seq.size()-1;
        codon::base hopping_base = this->seq[r_idx].pop(1);
        while (--r_idx) {
            hopping_base = this->seq[r_idx].squeeze_right(hopping_base);
        }
        this->seq[r_idx].insert_right(hopping_base);
    }
}

void codon::Seq::right_shift() {
    //removes right-most base in the first codon and start squeeze chain to the back propogating
    //bases until the final one --> if last codon is already full a new one will be generated,
    //increasing codon::Seq::seq.size() by one
    std::size_t idx{0};
    codon::base hopping_base = this->seq[idx].pop(this->seq[idx].get_bases_len());
    while (++idx) {
        hopping_base = this->seq[idx].squeeze_left(hopping_base);
    }
    if (this->seq.at(this->seq.size()-1).get_bases_len() < 3) {
        this->seq[idx].insert_left(hopping_base);
    } else {
        hopping_base = this->seq[idx].squeeze_left(hopping_base);
        this->seq.emplace_back(codon::Codon(hopping_base));
        //this should be fine assuming you have enough buffer left.
        //if not the vector resizes automatically but it should only happen once.
    }


}

void codon::Seq::insert_base(codon::base base, std::size_t insert_loc, int shift_loc=0) {
    //insert_loc = [0, 1, ..., this->seq.size()-1] index of seq where insertions should be taking place
    //shift_loc  = [0, 1, 2] where insertion will happen in codon:
    //  [0] base_1 [1] base_2 [2] base_3 -- [3] is not necessary; just use 0 in next Codon]
    //  any number above 2 will be treated as 2, squeezing out base 3.

    if (this->seq[insert_loc].get_bases_len() < 3) {
        //incase insert_loc is already an incomplete codon
        switch(shift_loc) {
            case 0: this->seq[insert_loc].insert_left(base); break;
            case 1: if (this->seq[insert_loc].get_bases_len() == 2) {
                        codon::base temp = this->seq[insert_loc].pop(2);
                        this->seq[insert_loc].insert_right(base);
                        this->seq[insert_loc].insert_right(temp);
                    }
                    else this->seq[insert_loc].insert_right(base); 
                    break;
            case 2: this->seq[insert_loc].insert_right(base); break;

        }
        return;
    }

    codon::base hopping_base;

    if (shift_loc == 0) {
        codon::base hopping_base = this->seq[insert_loc].squeeze_left(base);
    } else if (shift_loc == 1) {
        codon::base hopping_base = this->seq[insert_loc].pop(3);
        codon::base temp = this->seq[insert_loc].pop(2);
        this->seq[insert_loc].insert_right(base);
        this->seq[insert_loc].insert_right(temp);
    } else {
        codon::base hopping_base = this->seq[insert_loc].pop(3);
        this->seq[insert_loc].insert_right(base);
    }
    ++insert_loc;

    while(insert_loc != this->seq.size()-1) {
        if (this->seq[insert_loc].get_bases_len() == 3)
            this->seq[insert_loc].squeeze_left(hopping_base);
        else {
            this->seq[insert_loc].insert_left(hopping_base);
            return;
            // no need to propogate anymore
        }
        ++insert_loc;
    }
    //can reach this point if final Codon is already full and we have to make a new codon
    //can lead to resizing but effect is minimal because of existing buffer
    this->seq.emplace_back(codon::Codon(hopping_base));
}


void codon::Seq::insert_codon(codon::Codon codon_insert, std::size_t insert_loc, int shift_loc) {
    //insert_loc = [0, 1, ..., this->seq.size()-1] index of seq where insertions should be taking place
    //shift_loc  = [0, 1, 2] where insertion will happen in codon:
    //  [0] base_1 [1] base_2 [2] base_3 --> [3] is not necessary; just use 0 in next Codon]
    //  any number above 2 will be treated as 2, squeezing out base 3.

    int size_main = this->seq[insert_loc].get_bases_len();
    int size_other= codon_insert.get_bases_len();
    //edge case: insert_loc and insert can fit in the already existing codon
    if (size_main + size_other < 3){
        while (size_other--) {
            //popping and insert from reverse so that insert_loc can be reused for all insertions.
            this->insert_base(codon_insert.pop(codon_insert.get_bases_len()), insert_loc, shift_loc);
            return;
        }
    }
    //make space and new buffer if completely full
    if ((this->seq.size() + 1) < this->seq.capacity()) {
        this->seq.reserve((this->seq.size() + 1)*1.2);
        PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
    }

    //STEP 1 REARRANGE AND COMBINE - determine how much is right of shift_loc
    int amount_switch = this->seq[insert_loc].get_bases_len() - shift_loc;
    codon::Codon temp_reversed = Codon("VOID");
    while (amount_switch--) {
        temp_reversed.insert_right(this->seq[insert_loc].pop(0));
    }
    //fill up original insert_loc codon
    while (this->seq[insert_loc].get_bases_len() < 3) {
        if (codon_insert.get_bases_len() > 0) 
            this->seq[insert_loc].insert_right(codon_insert.pop(1));
        if (temp_reversed.get_bases_len() > 0) 
            codon_insert.insert_right(temp_reversed.pop(0));
    }
    while(temp_reversed.get_bases_len() > 0) 
        codon_insert.insert_right(temp_reversed.pop(0));

    // STEP 2 PUSH THAT INSERT IN (at this point we should be complete with the original insert_loc)
    switch (codon_insert.get_bases_len()) {
        case 1:
            this->insert_base(codon_insert.get_base_at(1), insert_loc + 1); 
            break;
        case 2: {
            //There is definitely a better way to do this but this will do for now
            std::vector<codon::Codon>::iterator it_seq = this->seq.begin()+insert_loc+1;
            this->seq.insert(it_seq, codon_insert);

            std::size_t rev_idx{this->seq.size()-1};
            codon::base hopping_base = this->seq[rev_idx].pop(1);

            while (--rev_idx > insert_loc+1) {
               hopping_base = this->seq[rev_idx].squeeze_right(hopping_base);
            }
            this->seq[rev_idx].insert_right(hopping_base); 
            return; 
        }
        case 3:
            std::vector<codon::Codon>::iterator it_seq = this->seq.begin()+insert_loc+1;
            this->seq.insert(it_seq, codon_insert);
            return;
    }
}


void codon::Seq::insert_seq(codon::Seq other, std::size_t insert_loc, int shift_loc=0) {
    //check if capacity is compromised
    if ((this->seq.size() + other.seq.size()) < this->seq.capacity()) {
        this->seq.reserve((this->seq.size() + other.seq.size())*1.2);
        PLOGD << "RESERVING MORE MEMORY FOR SEQUENCE";
    };
    int bases_at_insert = this->seq[insert_loc].get_bases_len();

    //REST NOT YET WORKED OUT YET
}


