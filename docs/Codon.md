
# codon documentation: Codon
## Introduction

Codon is the base class used in the codon library.
Each object contains a single byte (std::uint8_t) of encoded nucleotides and it supports 1 - 3 nucleotides.

Each pair of bits encodes for one location usable in the encoding process and can be visualised as follows:
> 00 00 00 00 = bits
> 00 11 22 33 = available loci

A codon is read left to right and at first the start ist evaluated by a marker => 01 or 10, which signifies the start of nucleotides inside the codon. Conversely, any 00 or 11 before combinations before a marker are ignored .
> __01__ 00 00 00 = Marker at LOC0
> 00 __10__ 00 00 = Marker at LOC1
> 11 11 __01__ 00 = Marker at LOC2
> ~~00 00 00 __01__ = Marker at LOC3~~ // invalid (would leave no space for any nucleotides)

Once the marker has been determined the remaining bit pairs to the right of the markers can be read as nucleotides A, G, C and T.

> A = 00
> G = 01
> C = 10
> T = 11

A few examples:
> AGT = 01 00 01 11
> TA = 11 10 11 00
> CAG = 01 01 00 01
> A = 11 11 10 00
> A = 00 00 01 00


## Inherent complement

You might have noticed that there exist a certain redundancy in the last example and in the design of the markers. The reason for this is that the marker not only signal the start of the nucleotides; it also contains information of about the orientation of this specific codon.

> 01 Marker = 5' => 3' (or N => C terminus for proteins)\
> 10 Marker = 3' => 5' (or C => N terminus for proteins)

A major advantage of this is that to get the complementary part of a nucleotide strand, only a single bitwise operation (~) is needed in order to flip the bits and the codon.

Consequently, this means that our previous examples can also be read like this:
> __01__ 00 01 11 = 5'-AGT-3'
> 11 __10__ 11 00 = 3'-TA-5'
> __01__ 01 00 01 = 5'-CAG-3'
> 11 11 __10__ 00 = 3'-A-5'
> 00 00 __01__ 00 = 5'-A-3'

The bits left of the marker are ignored but usually follow the same rule with 00 hinting at a 5'=>3' orientation and 11 vice versa. This, however, is not always given and if by some reason a codon like the following is generated it might catch you off guard:
> 11 __01__ 11 01 = 5'-TG-3'

### bp location / shift
In order to understand a couple of other functionalities it is important to know how encoded nucleotides are numbered, which is left to right, starting at 1.


| Codon     | binary representation | location / shift    |
| :--------:|:---------------------:|--------------------:|
| 5'-AGT-3' | 01 __00 01 11__       | A = 1, G = 2, T = 3 |
| 3'-TA-5'  | 11 10 __11 00__       | T = 1, A = 2        |
| 5'-CAG-3' | 01 __01 00 01__       | C = 1, A = 2, G = 3 |
| 3'-A-5'   | 11 11 10 __00__       | A = 1               |
| 5'-A-3'   | 00 00 01 __00__       | A = 1               |

Why not 0? There is no base pair 0 ... Empty is empt-

## VOID

That being said, the codon class also supports an object of 0 length, meaning an empty Codon with no nucleotides.
>5'-VOID-3' = 00 00 00 00
>3'-VOID-5' = 11 11 11 11

This object usually appears during a bug or a wrong implementation unless specifically created by the user (E.g. you want to declare and initialise a Codon object ahead of time but do not know what will end up being placed inside).

Most methods will throw if they are being passed a VOID but you can still insert into a VOID effectively turning it into a normal Codon.

VOIDs have a few pecularities in combination with more advanced features that you should watch out for:
* Codon::get_base_lenght = 0 but in a sequence, Seq::get_seq_len() will account for it when calculating the amount of stored objects.
* A VOID sequence (a Seq Object with only one or more VOIDs) behaves abnormally:
    * Returned locators are invalid:
        * get_first_loc() return an invalid locator past the sequence similar to .end() with iterators.
        * get_last_loc() returns {0, 1} (see more about locators in the respective chapter).

Don't use VOID codons unless absolutely necessary.

## Encoded characters
Since codons are stored in one byte it also sometimes makes sense to use an encoded character to represent the object in ascii compatible format. If cast directly this will most often not result in any readable output. The codon lib provides methods for turning the codons into printable characters and back again. This is also what is stored inside a .codon file instead of the nucleotide sequence to save space.

[//]: # (TODO: Add more information about encoded characters)

## Constructor
Codons can be constructed in various ways, see below

| Construct from    | Use case                                          |
| ----------------  |:--------------------------------------------------|
| std::string_view  | The standard use case                             |
| char              | Reserved for encoded characters                   |
| codon::base       | Creating a single base codon from an enum         |
| const Codon* const| will copy an existing Codon passed via pointer    |
| const Codon&      | will copy an existing Codon passed via reference  |
| Codon&&           | will move an existing Codon                       |

> !! The constructer from characters is reserved for encoded symbols, use the base enum for single base codons:
```c++
  codon::Codon not_what_you_think1{'A'};    // AAG
  codon::Codon not_what_you_think2{'G'};    // AGA
  codon::Codon not_what_you_think3{'C'};    // ACA
  codon::Codon not_what_you_think4{'T'};    // GGG

  codon::Codon doing_it_right1{codon::base::A}; // A
  codon::Codon doing_it_right2{codon::base::G}; // G
  codon::Codon doing_it_right3{codon::base::C}; // C
  codon::Codon doing_it_right4{codon::base::T}; // T
```

## Enumerators

Codon provides a few enumerators with the first one being the most commonly used.

#### codon::base : std::uint8_t
> codon::base::A = 0x00
> codon::base::G = 0x01
> codon::base::C = 0x10
> codon::base::T = 0x11

> codon::marker::n_strand_VOID = 0b00 00 00 00
> codon::marker::n_strand_1bp = 0b00 00 01 00
> codon::marker::n_strand_2bp = 0b00 01 00 00
> codon::marker::n_strand_3bp = 0b01 00 00 00
> codon::marker::c_strand_3bp = 0b10 00 00 00
> codon::marker::c_strand_2bp = 0b11 10 00 00
> codon::marker::c_strand_1bp = 0b11 11 10 00
> codon::marker::c_strand_VOID = 0b11 11 11 11

#### codon::mask : unsigned int
> codon::mask::base_1 = 00 00 00 11
> codon::mask::base_2 = 00 00 11 00
> codon::mask::mark_1 = 00 00 11 00 (= alias for base_2)
> codon::mask::r_half = 00 00 11 11
> codon::mask::base_3 = 00 11 00 00
> codon::mask::mark_2 = 00 11 00 00 (= alias for base_3)
> codon::mask::base_3 = 11 00 00 00
> codon::mask::l_half = 11 11 00 00

## Methods

### boolean verifiers
> is_full() -> bool;
> is_empty() -> bool;
> is_complement() -> bool;

```c++
codon::Codon empty_codon{"VOID"}; // Didn't you just- "shhhh sweet child"
codon::Codon semi_codon{"TA"};
codon::Codon full_codon{"GAC"};

empty_codon.is_full() // false
empty_codon.is_empty() // true

semi_codon.is_full() // false
semi_codon.is_empty() // false

full_codon.is_full() // true
full_codon.is_empty() // false

full_codon.is_complement() // false
full_codon.flip()
full_codon.is_complement() // true
```
### getters

> get_bases_len() -> int;
> get_bases_bin() -> std::bitset<8>;
> get_bases_str() -> std::string;
> get_bases_int() -> int;
> get_bases_encoded() -> char;
> get_base_at() -> codon::base;

```c++
codon::Codon start_codon{"AGT"};
start_codon.get_bases_len() // 3
start_codon.get_bases_bin() // 01000111
start_codon.get_bases_str() // AGT
start_codon.get_bases_int() // 71
start_codon.get_bases_encoded() // ;
start_codon.get_base_at(2) // codon::base::G
```

### inserters

> insert_right(codon::base);
> insert_left(codon::base);
> squeeze_right(codon::base) -> codon::base;
> insert_left(codon::base) -> codon::base;

```
codon::Codon start_codon{"AGT"};
codon::Codon semi_codon{"CA"};
start_codon.insert_right(codon::base::T); // ERROR ! Codon is full
codon::Codon dropped_base = start_codon.squeeze_right(codon::base::T);
// AGT + T => A + GTT
// dropped_base = codon::base::A
// start_codon = GTT
start_codon.squeeze_left(dropped_base); // A + GTT => AGT + T
semi_codon.insert_left(dropped_base); // A + CA => ACA

// insert ... only works if codon is not full
// squeeze ... only works if codon is full - drops/returns base on opposite end

```

### funky moves

> flip() -> codon::Codon;
> reverse() -> codon::Codon;
> flip_inplace();
> reverse_inplace();

```c++

codon::Codon start_codon{"AGT"};
codon::Codon flipped_codon = start_codon.flip(); // AGT -> TCA
codon::Codon reversed_codon = start_codon.reverse(); // AGT -> TGA

flipped_codon.flip_inplace(); // TCA -> AGT
reversed_codon.reverse_inplace(); // TGA -> AGT
```

## Functions

> base_to_str(const codon::base&) -> char;

```c++
codon::base_to_str(codon::base::A) // 'A'
// not a method!
```

## Overloaded operators
> = copy assignment
> = move assignment
> == equality
> != inequality

```c++
codon::Codon start_codon{"AGT"};
codon::Codon copy_start = start_codon;
codon::Codon flipped_start = std::move(start_codon.flip());

start_codon == flipped_start // false
start_codon == copy_start    // true
start_codon != flipped_start // true
```
