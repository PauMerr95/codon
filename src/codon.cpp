#include "codon.h"
#include <cstdint>
#include <string>
#include <bitset>
#include <plog/Log.h>

constexpr std::uint8_t MARKER_LOC3_5 = static_cast<uint8_t>(0b01000000);
constexpr std::uint8_t MARKER_LOC3_3 = static_cast<uint8_t>(0b10000000);
constexpr std::uint8_t MARKER_LOC3   = static_cast<uint8_t>(0b11000000);

constexpr std::uint8_t MARKER_LOC2_5 = static_cast<uint8_t>(0b00010000);
constexpr std::uint8_t MARKER_LOC2_3 = static_cast<uint8_t>(0b00100000);
constexpr std::uint8_t MARKER_LOC2   = static_cast<uint8_t>(0b00110000);

constexpr std::uint8_t MARKER_LOC1_5 = static_cast<uint8_t>(0b00000100);
constexpr std::uint8_t MARKER_LOC1_3 = static_cast<uint8_t>(0b00001000);
constexpr std::uint8_t MARKER_LOC1   = static_cast<uint8_t>(0b00001100);

/* For extraction purposes:
 * MARKER_LOC2 == base 1 in triplet
 * MARKER_LOC1 == base 2 in triplet
 * enum base T == base 3 in triplet
 *
 * base 1 is left most in the bit representation
 */

constexpr std::uint8_t DEL_LEFT_SIDE = static_cast<uint8_t>(0b00001111);

constexpr std::uint8_t VOID_5        = static_cast<uint8_t>(0b00000000);
constexpr std::uint8_t SWITCH_5      = static_cast<uint8_t>(0b11111111);


codon::Codon::Codon(const std::string& bases_str) {
    /* This function builds the bases from string using a 16bit generator.
     * This is probably not necessary but it was one of the things I added during debugging and 
     * it's only a temporary object.
     */
    if (bases_str == "VOID") {
        this->bases = VOID_5;
    }
    else if (bases_str == "SWITCH") {
        this->bases = SWITCH_5;
    }
    else {
        std::uint16_t generator {0};
        generator |= G;
        for (int i=0; i < bases_str.length(); i++) {
            switch (bases_str[i]) {
                case 'A': generator = generator << 2 | A; break;
                case 'G': generator = generator << 2 | G; break;
                case 'C': generator = generator << 2 | C; break;
                case 'T': generator = generator << 2 | T; break;
            }
        }
        this->bases = static_cast<uint8_t>(generator);
    }
};

codon::Codon::Codon(base base){
    this->bases = MARKER_LOC1_5 | base;
}

codon::Codon::~Codon() {
}


std::bitset<8> codon::Codon::get_bases_bin() const{
	return std::bitset<8>(this->bases);
}

int codon::Codon::get_bases_int() const {
	return static_cast<int>(this->bases);
}



/* This function returns the length of the inserted bases.
 * For void and switch codons it returns 0;
 */
std::size_t codon::Codon::get_bases_len() const {

    if (this->bases == VOID_5 || this->bases == SWITCH_5) return 0;

	int marker_3 = static_cast<std::uint8_t>(this->bases & MARKER_LOC3);
	int marker_2 = static_cast<std::uint8_t>(this->bases & MARKER_LOC2);

	if 	    (marker_3 == MARKER_LOC3_5 || marker_3 == MARKER_LOC3_3) return 3;
	else if (marker_2 == MARKER_LOC2_5 || marker_2 == MARKER_LOC2_3) return 2;
	else    return 1;
}



/* This function returns the bases as string and can be used for displaying.
 * It is reliant on get_bases_len() to calculate the length of the codon.
 * For void and switch codons the function aborts early return an empty string.
 */
std::string codon::Codon::get_bases_str() const{
    std::size_t len = this->get_bases_len();

    if (len == 0)
        return (this->bases == VOID_5) ? "VOID" : "SWITCH";

    std::string codon_str{};
    codon_str.resize(len);
    int idx {0};

    while (len) {
        std::uint8_t extracted_bits_at_len 
            = static_cast<std::uint8_t>(T) << (len-1)*2 & this->bases;

        //Shift into position 1 and put into switch statement:
        switch (static_cast<std::uint8_t>(extracted_bits_at_len >> (len-1)*2)) {
            case A: codon_str[idx] = 'A'; break;
            case G: codon_str[idx] = 'G'; break;
            case C: codon_str[idx] = 'C'; break;
            case T: codon_str[idx] = 'T'; break;
            default: 
                PLOGF << "Fatal error: Extracted bits could not be evaluated to A, G, C, T";
        }
        --len;
        ++idx;
    }

    return codon_str;
}



void codon::Codon::cast_to_switch() {
    if      (this->bases == VOID_5)     this->bases = SWITCH_5;
    else if (this->bases == SWITCH_5)   this->bases = VOID_5;
    else    PLOGF << "Fatal Error: Trying to cast encoding codon to switch";
}



void codon::Codon::insert_right(base base) {
    /* left -> new position 3
     * no length check necessary because that has to be done before calling the fn
     */
    this->bases = this->bases << 2;
    this->bases |= base;
}



void codon::Codon::insert_left(base base) {
    /* left -> new position 1 */
    if (this->get_bases_len() == 2) {
        this->bases &= DEL_LEFT_SIDE;
        this->bases |= (base << 4) | MARKER_LOC3_5;
    }
    else {
        this->bases &= static_cast<uint8_t>(T);
        this->bases |= (base << 2) | MARKER_LOC2_5;
    }
}

// Pushes base into base 3 returning previous base 1
codon::base codon::Codon::squeeze_right(base new_base) {
    enum base dropped_base = static_cast<enum base>((MARKER_LOC2 & this->bases) >> 4);
    //MARKER_LOC2 == BASE 1 for triplet, shifted by 4times so it can be converted to base
    this->bases <<= 2;
    this->bases |= new_base;
    this->bases &= ~(MARKER_LOC3);
    this->bases |= MARKER_LOC3_5;

    return dropped_base;
}

// Pushes base into base 1 returning previous base 3
codon::base codon::Codon::squeeze_left(base new_base) {
    enum base dropped_base = static_cast<enum base>((this->bases & T));
    this->bases >>= 2;
    this->bases &= DEL_LEFT_SIDE;
    this->bases |= static_cast<uint8_t>(new_base << 4) | MARKER_LOC3_5;

    return dropped_base;
}
