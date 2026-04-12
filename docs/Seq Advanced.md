
# codon documentation: Seq Advanced

## Good Job

Good Job making past locators. Knowing how they work is essential for most of the Seq methods and operations.

## Getters

### Get positionals:
```c++
codon::Seq sequence{"AGTTAACAGGATACGA"};
// AGT TAA CAG GAT ACG A
// ^-- --- --- --- --- ^
// -0- -1- -2- -3- -4- -5-
// 123 --- --- --- --- 1

codon::locator start_loc = sequence.get_first_loc();
// {0, 1}
codon::locator end_loc = sequence.get_last_loc();
// {5, 1}

std::size_t start_idx = sequence.get_first_idx(); // 0
std::size_t final_idx = sequence.get_first_idx(); // 5
```

### Get std::string

```c++
codon::Seq sequence{"AGTTAACAGGATACGA"};
// AGT TAA CAG GAT ACG A

// Two versions
std::string seq_string = sequence.get_seq_str(); // AGTTAACAGGATACGA
std::string seq_string_sep = sequence.get_seq_strsep(); // AGT TAA CAG GAT ACG A

These two methods are usually called without arguments but they both accept an output format:

sequence.get_seq_str(codon::OutputFormat::as_DNA);
// AGTTAACAGGATACGA

sequence.get_seq_str(codon::OutputFormat::as_RNA);
// AGUUAACAGGAUACGA

sequence.get_seq_str(codon::OutputFormat::as_PROT_w);
// SerTerGlnAspThr (The last 'A' is discarded because only a full codon can decode for protein)

sequence.get_seq_str(codon::OutputFormat::as_PROT);
// SXQDT (The last 'A' is discarded because only a full codon can decode for protein)

sequence.get_seq_str(codon::OutputFormat::as_CDN);
// (The last 'A' is NOT discarded)

// get_seq_strsep() also allows you to specify the seperator:
sequence.get_seq_strsep(codon::OutputFormat::as_PROT_w, '-');
// Ser-Ter-Gln-Asp-Thr-
// Note that the separator is also attached to the end. This is on purpose because it allows for nothing special - I am just too lazy to change it at the moment.

But there is more. There is an alternative overload that allows you to specify a range, using - you guessed it - locators.
std::pair<codon::locator, codon::locator> segment = {
    sequence.get_first_loc(),
    sequence.get_last_loc() - 1 // lets remove that pesky trailing 'A'
};

sequence.get_seq_strsep(segment, codon::OutputFormat::as_RNA);
// AGU UAA CAG GAU ACG

```

### Get as a Seq or Codon

```c++
codon::Seq sequence{"AGTTAACAGGATACGA"};
// AGT TAA CAG GAT ACG A

// Get a Seq
codon::Seq first_two_codons{sequence.subseq({0, 1}, {1, 3})};
// or
codon::locator start(0, 1);
codon::locator start(1, 3);
codon::Seq first_two_codons{sequence.subseq(start, end)};

// Get a Codon
codon::Codon first_codon{sequence.get_codon_at(sequence.get_first_loc())};
// AGT
// The first parameter is a non optional locator.
// The second parameter is an optional size 1 - 3: Default = 3.
// The third parameter is an optional boolean for codon overflow into the next: Default = false.
codon::locator start_cdn(2, 2);

// Sequence: AGT TAA CAG GAT ACG A
// Locator : --- --- -^- --- --- -

codon::Codon get_two{sequence.get_codon_at(start_cdn, 2)};
// AG
codon::Codon get_three{sequence.get_codon_at(start_cdn, 3)};
// AG ?

// Sequence: AGT TAA CAG GAT ACG A
// Locator : --- --- -^- --- --- -
// Locator : --- --- -12 3-- --- - Third one requires overflow to be set to true

codon::Codon get_three{sequence.get_codon_at(start_cdn, 3, true)};
// AGG

```

### Getting a length

Take the following example ...

```c++
codon::Seq sequence{AGTAACAGGATACGA};
// Sequence: AGT AAC AGG ATA CGA
// Locators: --- --- --- --- ---

codon::Seq sequence_shifted{AGTAACAGGATACGA};
sequence_shifted.right_shift();
// Sequence: AG TAA CAG GAT ACG A
// Locators: -- --- --- --- --- -

```

These two sequences are the same. But they have a different structure and those a different length depending on how you view them. That is why codon offers three different ways to access the length.
> .get_seq_len();

> .get_seq_trulen("codons");

> .get_seq_trulen("bp");

```c++
std::Seq VOID_SEQ {codon::Codon("VOID")}; // 👻

// Lenght of the underlying vector
sequence.get_seq_len(); // 5
sequence_shifted.get_seq_len(); // 6
VOID_SEQ.get_seq_len(); // 1

// Length of codons
sequence.get_seq_trulen(); // 5
sequence_shifted.get_seq_trulen(); // 6
VOID_SEQ.get_seq_trulen(); // 0 - not counting VOID

// Length of bases
sequence.get_seq_trulen("bp"); // 15
sequence_shifted.get_seq_trulen("bp"); // 15
VOID_SEQ.get_seq_trulen("bp"); // 0 - not counting VOID

// Checking the lenght of codons equals to the default value "codons"
sequence.get_seq_trulen() == sequence.get_seq_trulen("codons");
```

## Sequence manipulators

In the previous example we already encountered .right_shift(), a shifter that moved the sequence by one base to the right.

### .right_shift() | .left_shift()

```c++
codon::Seq sequence{AGTAACAGGATACGA};
// Sequence: AGT AAC AGG ATA CGA
// Locators: --- --- --- --- ---

sequence.right_shift(); // always in_place
// Sequence: AG TAA CAG GAT ACG A
// Locators: -- --- --- --- --- -

sequence.left_shift(); // always in_place
// Sequence: AGT AAC AGG ATA CGA
// Locators: --- --- --- --- ---
```
---
**NOTE**

    Shifters also accept an optional "up to index". This is primarily used within the library to close gaps with little use outside because creating a gap conventially should not be possible. It is mentioned here for completeness.Example incomplete index 2:

```c++
Sequence: AGT AAC AG ATA CGA
Locators: --- --- -- --- ---
sequence.left_shift(2);\
Sequence: AGT AAC AGA TAC GA
```

    Both shifters will also throw if the "up to index" is already full (the only exception being the final codon). This also means the following example cannot be left shifted
```c++
> Sequence: AGT AAC AGA TAC GA\
You can achieve the same effect by right shifting twice\
> Sequence: A GTA ACA GAT ACG A
```
---

### .reverse() | .reverse_inplace()

Reverses the entire or a fragment of the sequence, either in place or by returning the modified sequence.

```c++
// Sequence: AGT AAC AG ATA CGA
// Locators: --- --- -- --- ---

sequence.reverse();
// Sequence: AGC ATA GGA CAA TGA
// Locators: ^-- --- --- --- --^

sequence.reverse({2, 2}, {4, 1});
// Sequence: AGT AAC ACA TAG GGA
// Locators: --- --- -^- --- ^--
```

### .flip() | .flip_inplace()

Flips the entire or a fragment of the sequence, either in place or by returning the modified sequence.

```c++
// Sequence: AGT AAC AG ATA CGA
// Locators: --- --- -- --- ---

sequence.flip();
// Sequence: TCG TAT CCT GTT ACT
// Locators: ^-- --- --- --- --^

sequence.flip({2, 2}, {4, 1});
// Sequence: AGT AAC ACT GTT AGA
// Locators: --- --- -^- --- ^--
```

## Removers

### pop_base(codon::locator) -> codon::base

```c++
// Sequence: AGT AAC AGA TAC GA
// Locators: --- --- -^- --- --

sequence.pop_base(codon::locator(2, 2)); // codon::base::G
// Sequence: AGT AAC AAT ACG A
// Gap is automatically filled
```

### pop_codon(codon::locator, int) -> codon::Codon

Additionally accepts a value for the size to be extracted. Default value = 3.

---
**NOTE**

    In case the size is larger than what is capable of being removed the maximum possible amount will be removed and returned.

---

```c++
// Sequence: AGT AAC AGA TAC GA
// Locators: --- --- -^- --- --

sequence.pop_codon(codon::locator(2, 2), 2); // returns codon::Codon("GA")
// Sequence: AGT AAC ATA CGA
// Locators: --- --- --- ---
// Gap is automatically filled
```

### pop_seq(codon::locator, std::size_t) -> codon::Seq
### pop_seq(codon::locator, codon::locator) -> codon::Seq

```c++
// Sequence: AGT AAC AGA TAC GA
// Locators: --- ^-- --- --- --
// Locators: --- 123 45- --- --

sequence.pop_codon(codon::locator(1, 1), 5);
// returns codon::Seq("AACAG") - automatically left aligns
// Sequence: AGT ATA CGA
// Locators: --- --- ---
// Gap is automatically filled
```

## Extenders

There are two types of extenders: inserters and an overloaded push_back

### insert_base(codon::base, codon::locator)
```c++
// Sequence: AGT AAC AGA TAC GA
// Locators: --- ^-- --- --- --

sequence.insert_base(codon::base::T, codon::locator(1, 1));
// Sequence: AGT TAA CAG ATA CGA
// Insert  : --- ^-- --- --- --

```
### insert_codon(codon::Codon, codon::locator)
```c++
// Sequence: AGT AAC AGA TAC GA
// Locators: --- --- -^- --- --

sequence.insert_base(codon::Codon("AGA"), codon::locator(2, 2));
// Sequence: AGT TAA CAG AAG ATA CGA
// Insert  : --- --- -^^ ^-- --- ---
```
### insert_seq(codon::seq, codon::locator)
```c++
// Sequence: AGT AAC AGA TAC GA
// Locators: --- --^ --- --- --

sequence.insert_base(codon::Seq("AGCAGGC"), codon::locator(1, 3));
// Sequence: AGT AAA GCA GGC CAG AAG ATA CGA
// Insert  : --- --^ ^^^ ^^^ --- --- --- ---
```

---
**NOTE**

    Only the shift of the right part after the insertion locus is affected by the insertion.

---
### push_back(codon::base)
```c++
// Sequence: AGT AAC AGA TAC GA
sequence.push_back(codon::base::T);
// Sequence: AGT AAC AGA TAC GAT
// Insert  : --- --- --- --- --^
```
### push_back(codon::Codon)
```c++
// Sequence: AGT AAC AGA TAC GA
sequence.push_back(codon::Codon("CAA"));
// Sequence: AGT AAC AGA TAC GAC AA
// Insert  : --- --- --- --- --^ ^^
```
### push_back(codon::Seq)
```c++
// Sequence: AGT AAC AGA TAC GA
sequence.push_back(codon::Seq("AGGCAA"));
// Sequence: AGT AAC AGA TAC GAA GGC AA
// Insert  : --- --- --- --- --^ ^^^ ^^
```

---
**NOTE**

    Depending on the shift of the locator the alignment of the extending item will be affected.

---
