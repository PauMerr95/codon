#include "codon.h"
#include "seq.h"
#include "testing.h"
#include <cassert>
#include <plog/Log.h>
#include <string>


int test::seq_test() {
    const std::string test_seq_1_str =
        "AGCTTGACGATGATCGATTTCGAACTGGCATGGGACAGTACTAGCATAGCATGCTAGCTGGATCGACT";
    const std::string test_seq_2_str =
        "CGAACTGGCATGGGACAGTACTAGCATAGCATGCTAGCTGGATGACTAGCTTGACGATGATCGATTT";
    const std::string test_seq_3_str =
        "CGAACTGGCATGGGACCTAGCTTGACGATGATCGATTTAGTACTACATAGCATGCTAGCTGGATCGA";
    const std::string test_seq_4_str =
        "GTACTAGCATAGCATGCTAGCTGGATCGACTAGCTTGACGATGATCGTCGAACTGGCATGGGACAT";

    codon::Seq test_seq_1 = codon::Seq(test_seq_1_str);
    assert(test_seq_1.get_encoding_str() == test_seq_1_str);

    codon::Seq test_seq_2 = codon::Seq(test_seq_2_str);
    assert(test_seq_2.get_encoding_str() == test_seq_2_str);

    codon::Seq test_seq_3 = codon::Seq(test_seq_3_str);
    assert(test_seq_3.get_encoding_str() == test_seq_3_str);

    codon::Seq test_seq_4 = codon::Seq(test_seq_4_str);
    assert(test_seq_4.get_encoding_str() == test_seq_4_str);

    return 0;
}


#ifdef SEQ_TEST

int main () {

    assert(test::codon_test() == 0);
    PLOGD << "Passed sequence test";
    return 0;
}

#endif /* ifdef  */
