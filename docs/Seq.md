
# codon documentation: Seq Basics and Locators

## Introduction

TLDR: Seq is a wrapper for a vector of strings.

Jokes aside, there is actually quite a big difference in comparison to normal strings. Because codons contain up to three bases each, the structuring / shift can and has to be considered.

Take the following string as an example:
```c++
> ATGTAACAGGAT
// There are three meaningful ways to represent this string in a Seq
> ATG TAA CAG GAT
>  AT GTA ACA GGA T
>   A TGT AAC AGG AT
```
Well, which one is correct? Technically they all are but it depends on your use case. Like in a real biological system, if you want to get a translation the Seq needs to be aligned properly. In this case, considering that ATG is a start codon only the first version makes sense. This, however can get increasingly complicated, especially after having manipulated your data.
(No worries, though, there are ways to find out quickly what alignment you want to store your sequence, that will be covered later - more on that in the chaper on Maps.)

## Creating a Seq

Not accounting for copy and move constructors there are three major ways to initiate a Seq.

### 1. Construct from string (basic or encoded)

```c++
// Full constructor:
// codon::Seq::Seq(std::string_view input, codon::format::AGCT)
// Takes two formats: "AGCT" and "encoded", default is AGCT
codon::Seq string_seq{"AGTTAACAGGAT"};
// AGT TAA CAG GAT

codon::Seq string_seq("&fe^h", codon::format::encoded);
// ACG TCG CGT TCC

// Will always create a version that is full on the 5' end.
// This constructor will not allow you to create a VOID sequence!
codon::Seq invalid{"VOID"}; // !! will throw
```

### 2. Construct empty with reserved size
```c++
// Full constructor:
// codon::Seq::Seq(std::size_t)

codon::Seq string_seq{10'000};

// {} but with reserved size of 10'000
// not the same as a VOID sequence
```


### 3. Construct from Codon or Seq
```c++
// Possible with copy or move.
codon::Codon not_a_compiler{"GCC"};

codon::Seq seq_codon_copy{not_a_compiler}; // copies the original
codon::Seq seq_seq_copy{seq_codon_copy}; // copies the previously created sequence

codon::Seq seq_codon_move(std::move(codon::Codon("CAC")));
// the temporary Codon is moved and the newly created Seq takes ownership

codon::Seq seq_seq_move(std::move(seq_codon_copy.flip()));
// the temporary flipped Seq is moved and the newly created Seq takes ownership

//Creating a Seq from a codon also allows for the following which is the only valid way to create a VOID Seq
codon::Seq doNotTryAtHome1(codon::Codon("VOID"));
codon::Seq doNotTryAtHome2(std::move(codon::Codon("VOID")));
```

## Locators

Locators are a way to reference a part of the sequence in a very simple form. Like iterator they are not automatically updated and can get invalidated. Essentially they are just aggregates containing an index for the codon in the sequence and a shift to point at the base.

As such only shifts between 1 and 3 are valid.
The index is a std::size_t or unsigned long long capable of storing very large positions.

A few examples:
```c++
> locator {3, 1}
> AGT CTA AAA ACC GAG GA
> ___ ___ ___ ^__ ___ __
> -0- -1- -2- -3- -4- -5- INDEX
> --- --- --- 123 --- --- SHIFT

> locator {0, 3}
> AGT CTA AAA ACC GAG GA
> __^ ___ ___ ___ ___ __
> -0- -1- -2- -3- -4- -5- INDEX
> 123 --- --- --- --- --- SHIFT

> locator {5, 2}
> AGT CTA AAA ACC GAG GA
> ___ ___ ___ ___ ___ _^
> -0- -1- -2- -3- -4- -5- INDEX
> --- --- --- --- --- 12- SHIFT

> locator {5, 3}
> AGT CTA AAA ACC GAG GA
> ___ ___ ___ ___ ___ __X INVALID LOCATOR
> -0- -1- -2- -3- -4- -5- INDEX
> --- --- --- --- --- 12- SHIFT
```
In case of incomplete codons the left most base is always 1, regardless of orientation.
Because a sequence can contain incomplete Codons it is not a given that only the last codon is potentially incomplete. Quite often the first one is as well to account for the global shift.
```c++
> locator {0, 3}
> AG- TCT AAA AAC CGA GGA
> __X ___ ___ ___ ___ ___ INVALID LOCATOR
> -0- -1- -2- -3- -4- -5- INDEX
> 123 --- --- --- --- --- SHIFT

> locator {5, 3}
> AG- TCT AAA AAC CGA GGA
> ___ ___ ___ ___ ___ __^ OKAY
> -0- -1- -2- -3- -4- -5- INDEX
> --- --- --- --- --- 123 SHIFT
```

// Constructing is simple:

```c++
codon::locator loc{2000, 2};

// Many functions return locators:
codon::Seq sample{"AGTCTAAAAACCGAGGA"};
codon::locator sample_begin = sample.get_first_loc;   // {0, 1}
codon::locator sample_end = sample.get_last_loc;      // {5, 2}

// Or expect locators:
sample.get_base_at(sample_begin); // codon::base::A

// Members are public and can be changed
sample_end.shift = 4;   // Do not do this!

// Locators contain a few methods:
sample_end.verify_shift();   // Throws because 4 is an invalid shift.
sample_end.shift = 3;
sample_end.verify_shift();   // Does not throw despite being an invalid locator

// to check the validity of a locator for a specific sequence you can use the following Seq method
sample.is_locator_valid(sample_begin);   // true
sample.is_locator_valid(sample_end);     // false

std::size_t distance = sample_begin.distance_to(sample_end);
// Return the distance (how many jumps separate the two locators).
// This is not a valid way to count bases because it assumes every codon inbetween is three bases long. Use codon::Seq::get_seq_trulen(codon::format::bp) instead

std::string locator_string = sample_end.to_str(); // "{5, 3}"

// Locators also support a lot of operators:
sample_begin == sample_end // false - requires both index and shift to be the same
sample_begin != sample_end // true - requires both index and shift to be the different

sample_begin < sample_end  // true - also accounts for shift differences
sample_begin > sample_end  // false - also accounts for shift differences

sample_begin <= sample_end // true - also accounts for shift differences
sample_begin >= sample_end // false - also accounts for shift differences

sample_begin + 5           // {0, 1} -> {1, 3} temporary expression of a locator that is 5 'jumps' ahead
sample_end - 5             // {5, 3} -> {4, 1} temporary expression of a locator that is 5 'jumps' behind
sample_end += 1            // {5, 3} will be now {6, 1}
sample_end -= 2            // {6, 1} will be now {5, 2}

```

__An important feature of locators is that they are end-inclusive, meaning that unlike iterators when passed to a method they will always include the final location. There is no traditional .end() that points past the final object!__
