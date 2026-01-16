#include "codon.h"
#include <cassert>
#include <bitset>
#include <plog/Log.h>


int main() {

	Codon triplet = Codon("TCA");
        assert(std::bitset<8>   (0b01000000) == triplet.get_bases_bin());
        assert(static_cast<int> (0b01000000) == triplet.get_bases_int());
        assert(triplet.get_bases_len() == 3);
        assert(triplet.get_bases_str() == std::string("TCA"));
    PLOGD << "Passed triplet check \n";

	Codon duplet = Codon("GC");
        assert(std::bitset<8>   (0b00010110) == duplet.get_bases_bin());
        assert(static_cast<int> (0b00010110) == duplet.get_bases_int());
        assert(duplet.get_bases_len() == 2);
        assert(duplet.get_bases_str() == std::string("GC"));
    PLOGD << "Passed duplex check \n";
	
	Codon singlet = Codon("T");
        assert(std::bitset<8>   (0b00000111) == singlet.get_bases_bin());
        assert(static_cast<int> (0b00000111) == singlet.get_bases_int());
        assert(singlet.get_bases_len() == 1);
        assert(singlet.get_bases_str() == std::string("T"));
    PLOGD << "Passed singlet check \n";

	return 0;
}
