#include "codon.h"
#include <cstdint>
#include <string>
#include <bitset>
#include <plog/Log.h>

#define MARKER_POS3_5   static_cast<uint8_t>(0b01000000)
#define MARKER_POS3_3   static_cast<uint8_t>(0b10000000)
#define MARKER_LOC3     static_cast<uint8_t>(0b11000000)

#define MARKER_POS2_5   static_cast<uint8_t>(0b00010000)
#define MARKER_POS2_3   static_cast<uint8_t>(0b00100000)
#define MARKER_LOC2     static_cast<uint8_t>(0b00110000)

#define MARKER_POS1_5   static_cast<uint8_t>(0b00000100)
#define MARKER_POS1_3   static_cast<uint8_t>(0b00001000)
#define MARKER_LOC1     static_cast<uint8_t>(0b00001100)

#define VOID_5          static_cast<uint8_t>(0b00000000) 
#define SWITCH_5        static_cast<uint8_t>(0b11111111) 


Codon::Codon(const std::string& bases_str) {
    /* This function builds the bases from string using a 16bit generator.
     * This is probably not necessary but it was one of the things I added during debugging and 
     * it's only a temporary object.
     */
    if (bases_str == "VOID") {
        this->bases = VOID_5;
        PLOGD << "Generated Codon " << bases_str << " as " << this->get_bases_bin();
    }
    else if (bases_str == "SWITCH") {
        this->bases = SWITCH_5;
        PLOGD << "Generated Codon " << bases_str << " as " << this->get_bases_bin();
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
        PLOGD << "Generated Codon " << bases_str << " as " << this->get_bases_bin();
    }
};

Codon::Codon(base base){
    this->bases = MARKER_POS1_5 | base;
    PLOGD << "Generated Codon from base: " << this->get_bases_bin();
}

Codon::~Codon() {
    PLOGD << "Destroying Codon at memory location " << &this->bases;
}


std::bitset<8> Codon::get_bases_bin() const{
	return std::bitset<8>(this->bases);
}

int Codon::get_bases_int() const {
	return static_cast<int>(this->bases);
}

//untested:
std::size_t Codon::get_bases_len() const {
    /* This function returns the length of the inserted bases.
     * For void and switch codons it returns 0;
     */

    PLOGD << "Running fn get_bases_len";

    if (this->bases == VOID_5 || this->bases == SWITCH_5) return 0;

	int marker_3 = static_cast<std::uint8_t>(this->bases & MARKER_LOC3);
	int marker_2 = static_cast<std::uint8_t>(this->bases & MARKER_LOC2);

	if 	    (marker_3 == MARKER_POS3_5 || marker_3 == MARKER_POS3_3) return 3;
	else if (marker_2 == MARKER_POS2_5 || marker_2 == MARKER_POS2_3) return 2;
	else    return 1;
}

//untested:
std::string Codon::get_bases_str() const{
    /* This function returns the bases as string and can be used for displaying.
     * It is reliant on get_bases_len() to calculate the length of the codon.
     * For void and switch codons the function aborts early return an empty string.
     */
    
    PLOGD << "Running fn get_bases_str";

    std::size_t len = this->get_bases_len();

    PLOGD << "Determined len to be " << len;

    if (len == 0) return ""; //For void and switch

    std::string codon_str{};
    codon_str.resize(len);
    int idx {0};

    while (len) {
        std::uint8_t extracted_bits_at_len 
            = static_cast<std::uint8_t>(A) << (len-1)*2 & this->bases;

        PLOGD << "Extracted bits at pos. " << idx << ": \t" 
            << std::bitset<8>(extracted_bits_at_len);

        //Shift into position 1 and put into switch statement:
        switch (static_cast<std::uint8_t>(extracted_bits_at_len >> (len-1)*2)) {
            case A: 
                PLOGD << "Extracted bits evaluated to 'A' and pushed to idx " << idx;
                codon_str[idx] = 'A'; break;
            case G: 
                PLOGD << "Extracted bits evaluated to 'G' and pushed to idx " << idx;
                codon_str[idx] = 'G'; break;
            case C: 
                PLOGD << "Extracted bits evaluated to 'C' and pushed to idx " << idx;
                codon_str[idx] = 'C'; break;
            case T: 
                PLOGD << "Extracted bits evaluated to 'T' and pushed to idx " << idx;
                codon_str[idx] = 'T'; break;
            default: 
                PLOGF << "Extracted bits could not be evaluated to A, G, C, T";
        }
        --len;
        ++idx;
    }

    PLOGD << "Determined string to be '" << codon_str << "' for " << this->get_bases_bin();

    return codon_str;
}

void Codon::cast_to_switch() {
    if      (this->bases == VOID_5)     this->bases = SWITCH_5;
    else if (this->bases == SWITCH_5)   this->bases = VOID_5;
    else    PLOGF << "Fatal Error: Trying to cast encoding codon to switch";
}

void Codon::insert_left(base base) {
    /* left -> new position 1
     * no length check necessary because that has to be done before calling the fn
     */
    this->bases = this->bases << 2;
    this->bases |= base;
    PLOGD << "Inserted base left at pos. 1, no dropped base";
}

void Codon::insert_right(base base) {
    /* left -> new position 2 or 3 */
    if (this->get_bases_len() == 2) {
        this->bases &= static_cast<uint8_t>(0b00001111);
        this->bases |= (base << 4) | MARKER_POS3_5;
    PLOGD << "Inserted base right at pos. 3, no dropped base";
    }
    else {
        this->bases &= static_cast<uint8_t>(0b00000011);
        this->bases |= (base << 2) | MARKER_POS2_5;
    PLOGD << "Inserted base left at pos. 2, no dropped base";
    }
}

base Codon::squeeze_left(base base) {
    /* left -> new pos. 1,
     * pos 3 is lost and will be return to be further passed down
     */
    this->bases = this->bases << 2;
    enum base dropped_base = static_cast<enum base>((this->bases & MARKER_LOC3) >> 6);
    this->bases &= static_cast<uint8_t>(0b00111111);
    this->bases |= MARKER_LOC3;

    PLOGD << "Squeezed base left at pos. 1 and dropped base " << dropped_base;
    return dropped_base;
}

base Codon::squeeze_right(base base) {
    /* left -> new pos. 3,
     * pos 1 is lost and will be return to be further passed down
     */
    enum base dropped_base = static_cast<enum base>((this->bases & T));
    this->bases &= static_cast<uint8_t>(0b00111111);
    this->bases |= MARKER_LOC3;

    PLOGD << "Squeezed base left at pos. 3 and dropped base " << dropped_base;
    return dropped_base;
}
