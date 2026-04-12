#include "transmute.h"

#include <unordered_map>

#include "codon.h"

std::unordered_map<codon::Codon, std::string> codon::translate_w{
    {codon::Codon("AAA"), "Lys"}, {codon::Codon("AAG"), "Lys"},
    {codon::Codon("AAC"), "Asn"}, {codon::Codon("AAT"), "Asn"},

    {codon::Codon("AGA"), "Arg"}, {codon::Codon("AGG"), "Arg"},
    {codon::Codon("AGC"), "Ser"}, {codon::Codon("AGT"), "Ser"},

    {codon::Codon("ACA"), "Thr"}, {codon::Codon("ACG"), "Thr"},
    {codon::Codon("ACC"), "Thr"}, {codon::Codon("ACT"), "Thr"},

    {codon::Codon("ATA"), "Ile"}, {codon::Codon("ATG"), "Met"},
    {codon::Codon("ATC"), "Ile"}, {codon::Codon("ATT"), "Ile"},

    {codon::Codon("GAA"), "Glu"}, {codon::Codon("GAG"), "Glu"},
    {codon::Codon("GAC"), "Asp"}, {codon::Codon("GAT"), "Asp"},

    {codon::Codon("GGA"), "Gly"}, {codon::Codon("GGG"), "Gly"},
    {codon::Codon("GGC"), "Gly"}, {codon::Codon("GGT"), "Gly"},

    {codon::Codon("GCA"), "Ala"}, {codon::Codon("GCG"), "Ala"},
    {codon::Codon("GCC"), "Ala"}, {codon::Codon("GCT"), "Ala"},

    {codon::Codon("GTA"), "Val"}, {codon::Codon("GTG"), "Val"},
    {codon::Codon("GTC"), "Val"}, {codon::Codon("GTT"), "Val"},

    {codon::Codon("CAA"), "Gln"}, {codon::Codon("CAG"), "Gln"},
    {codon::Codon("CAC"), "His"}, {codon::Codon("CAT"), "His"},

    {codon::Codon("CGA"), "Arg"}, {codon::Codon("CGG"), "Arg"},
    {codon::Codon("CGC"), "Arg"}, {codon::Codon("CGT"), "Arg"},

    {codon::Codon("CCA"), "Pro"}, {codon::Codon("CCG"), "Pro"},
    {codon::Codon("CCC"), "Pro"}, {codon::Codon("CCT"), "Pro"},

    {codon::Codon("CTA"), "Leu"}, {codon::Codon("CTG"), "Leu"},
    {codon::Codon("CTC"), "Leu"}, {codon::Codon("CTT"), "Leu"},

    {codon::Codon("TAA"), "Ter"}, {codon::Codon("TAG"), "Ter"},
    {codon::Codon("TAC"), "Tyr"}, {codon::Codon("TAT"), "Tyr"},

    {codon::Codon("TGA"), "Ter"}, {codon::Codon("TGG"), "Trp"},
    {codon::Codon("TGC"), "Cys"}, {codon::Codon("TGT"), "Cys"},

    {codon::Codon("TCA"), "Ser"}, {codon::Codon("TCG"), "Ser"},
    {codon::Codon("TCC"), "Ser"}, {codon::Codon("TCT"), "Ser"},

    {codon::Codon("TTA"), "Leu"}, {codon::Codon("TTG"), "Leu"},
    {codon::Codon("TTC"), "Phe"}, {codon::Codon("TTT"), "Phe"},
};

std::unordered_map<codon::Codon, char> codon::translate{
    {codon::Codon("AAA"), 'K'}, {codon::Codon("AAG"), 'K'},
    {codon::Codon("AAC"), 'N'}, {codon::Codon("AAT"), 'N'},

    {codon::Codon("AGA"), 'R'}, {codon::Codon("AGG"), 'R'},
    {codon::Codon("AGC"), 'S'}, {codon::Codon("AGT"), 'S'},

    {codon::Codon("ACA"), 'T'}, {codon::Codon("ACG"), 'T'},
    {codon::Codon("ACC"), 'T'}, {codon::Codon("ACT"), 'T'},

    {codon::Codon("ATA"), 'I'}, {codon::Codon("ATG"), 'M'},
    {codon::Codon("ATC"), 'I'}, {codon::Codon("ATT"), 'I'},

    {codon::Codon("GAA"), 'E'}, {codon::Codon("GAG"), 'E'},
    {codon::Codon("GAC"), 'D'}, {codon::Codon("GAT"), 'D'},

    {codon::Codon("GGA"), 'G'}, {codon::Codon("GGG"), 'G'},
    {codon::Codon("GGC"), 'G'}, {codon::Codon("GGT"), 'G'},

    {codon::Codon("GCA"), 'A'}, {codon::Codon("GCG"), 'A'},
    {codon::Codon("GCC"), 'A'}, {codon::Codon("GCT"), 'A'},

    {codon::Codon("GTA"), 'V'}, {codon::Codon("GTG"), 'V'},
    {codon::Codon("GTC"), 'V'}, {codon::Codon("GTT"), 'V'},

    {codon::Codon("CAA"), 'Q'}, {codon::Codon("CAG"), 'Q'},
    {codon::Codon("CAC"), 'H'}, {codon::Codon("CAT"), 'H'},

    {codon::Codon("CGA"), 'R'}, {codon::Codon("CGG"), 'R'},
    {codon::Codon("CGC"), 'R'}, {codon::Codon("CGT"), 'R'},

    {codon::Codon("CCA"), 'P'}, {codon::Codon("CCG"), 'P'},
    {codon::Codon("CCC"), 'P'}, {codon::Codon("CCT"), 'P'},

    {codon::Codon("CTA"), 'L'}, {codon::Codon("CTG"), 'L'},
    {codon::Codon("CTC"), 'L'}, {codon::Codon("CTT"), 'L'},

    {codon::Codon("TAA"), 'X'}, {codon::Codon("TAG"), 'X'},
    {codon::Codon("TAC"), 'Y'}, {codon::Codon("TAT"), 'Y'},

    {codon::Codon("TGA"), 'X'}, {codon::Codon("TGG"), 'W'},
    {codon::Codon("TGC"), 'C'}, {codon::Codon("TGT"), 'C'},

    {codon::Codon("TCA"), 'S'}, {codon::Codon("TCG"), 'S'},
    {codon::Codon("TCC"), 'S'}, {codon::Codon("TCT"), 'S'},

    {codon::Codon("TTA"), 'L'}, {codon::Codon("TTG"), 'L'},
    {codon::Codon("TTC"), 'F'}, {codon::Codon("TTT"), 'F'},
};

std::unordered_map<codon::Codon, std::string> codon::transcribe{
    {codon::Codon("AAA"), "AAA"}, {codon::Codon("AAG"), "AAG"},
    {codon::Codon("AAC"), "AAC"}, {codon::Codon("AAT"), "AAU"},

    {codon::Codon("AGA"), "AGA"}, {codon::Codon("AGG"), "AGG"},
    {codon::Codon("AGC"), "AGC"}, {codon::Codon("AGT"), "AGU"},

    {codon::Codon("ACA"), "ACA"}, {codon::Codon("ACG"), "ACG"},
    {codon::Codon("ACC"), "ACC"}, {codon::Codon("ACT"), "ACU"},

    {codon::Codon("ATA"), "AUA"}, {codon::Codon("ATG"), "AUG"},
    {codon::Codon("ATC"), "AUC"}, {codon::Codon("ATT"), "AUU"},

    {codon::Codon("GAA"), "GAA"}, {codon::Codon("GAG"), "GAG"},
    {codon::Codon("GAC"), "GAC"}, {codon::Codon("GAT"), "GAU"},

    {codon::Codon("GGA"), "GGA"}, {codon::Codon("GGG"), "GGG"},
    {codon::Codon("GGC"), "GGC"}, {codon::Codon("GGT"), "GGU"},

    {codon::Codon("GCA"), "GCA"}, {codon::Codon("GCG"), "GCG"},
    {codon::Codon("GCC"), "GCC"}, {codon::Codon("GCT"), "GCU"},

    {codon::Codon("GTA"), "GUA"}, {codon::Codon("GTG"), "GUG"},
    {codon::Codon("GTC"), "GUC"}, {codon::Codon("GTT"), "GUU"},

    {codon::Codon("CAA"), "CAA"}, {codon::Codon("CAG"), "CAG"},
    {codon::Codon("CAC"), "CAC"}, {codon::Codon("CAT"), "CAU"},

    {codon::Codon("CGA"), "CGA"}, {codon::Codon("CGG"), "CGG"},
    {codon::Codon("CGC"), "CGC"}, {codon::Codon("CGT"), "CGU"},

    {codon::Codon("CCA"), "CCA"}, {codon::Codon("CCG"), "CCG"},
    {codon::Codon("CCC"), "CCC"}, {codon::Codon("CCT"), "CCU"},

    {codon::Codon("CTA"), "CUA"}, {codon::Codon("CTG"), "CUG"},
    {codon::Codon("CTC"), "CUC"}, {codon::Codon("CTT"), "CUU"},

    {codon::Codon("TAA"), "UAA"}, {codon::Codon("TAG"), "UAG"},
    {codon::Codon("TAC"), "UAC"}, {codon::Codon("TAT"), "UAU"},

    {codon::Codon("TGA"), "UGA"}, {codon::Codon("TGG"), "UGG"},
    {codon::Codon("TGC"), "UGC"}, {codon::Codon("TGT"), "UGU"},

    {codon::Codon("TCA"), "UCA"}, {codon::Codon("TCG"), "UCG"},
    {codon::Codon("TCC"), "UCC"}, {codon::Codon("TCT"), "UCU"},

    {codon::Codon("TTA"), "UUA"}, {codon::Codon("TTG"), "UUG"},
    {codon::Codon("TTC"), "UUC"}, {codon::Codon("TTT"), "UUU"},
};
