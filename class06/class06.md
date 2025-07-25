# Class06: function
Kate Zhou (PID:A17373286)

- [1.Add function](#1add-function)
- [2. Generate DNA Function](#2-generate-dna-function)
- [3. Generate Protein Function](#3-generate-protein-function)

## 1.Add function

New function - name (function name) - input arguments (there can be
loads of these separated by a comma) - the body (the R code that does
the work)

``` r
add <- function(x, y=10, z = 100) {
  x+y+z
}
```

I can just use this function like any functions long as R knows about it
(i.e. run the code chunk)

``` r
add(1, 100)
```

    [1] 201

``` r
add(c(1, 2, 3, 4), 100)
```

    [1] 201 202 203 204

``` r
add(1)
```

    [1] 111

Functions can have “required” input arguments and “optional” input
arguments. The optional arguments are defined with an equals default
value (`y=100`) in the function definition

## 2. Generate DNA Function

> Q. Write a function to return a DNA sequence of a user specified
> length? Call it `generate_dna()`

``` r
students <- c("jeff", "jeremy", "peter")

sample(students, size=5, replace=TRUE)
```

    [1] "jeff"  "peter" "peter" "jeff"  "peter"

Now work with `bases` rather than `students`

``` r
bases <- c("A", "C", "G", "T")

sample(bases, size=10, replace=TRUE)
```

     [1] "T" "T" "T" "C" "T" "T" "A" "T" "G" "C"

Now I have a working “snipt” of

``` r
generate_dna <-function(size=5){
  bases <- c("A", "C", "G", "T")
  sample(bases, size=size, replace=TRUE)
}
```

``` r
generate_dna()
```

    [1] "G" "C" "G" "T" "G"

I want the ability to return a sequence like “AGCACCTG” (i.e a one
element vector where the bases are all together)

``` r
generate_dna <-function(size=5, together=TRUE) {
  bases <- c("A", "C", "G", "T")
  sequence <- sample(bases, size=size, replace=TRUE)
  if (together) {
    sequence <- paste(sequence, collapse = "")
  }
  return(sequence)
}
```

``` r
generate_dna()
```

    [1] "ACGGA"

## 3. Generate Protein Function

> Q. Write a protein sequence generating function that will return
> sequence of a user specified length?

We can get the set of 20 natural amino-acids from the **bio3d** package

``` r
library("bio3d")
bio3d::aa.table$aa1[1:20]
```

     [1] "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M" "F" "P" "S" "T" "W" "Y"
    [20] "V"

``` r
generate_protein <-function(size=6, together=TRUE) {
  bases <- bio3d::aa.table$aa1[1:20]
  
  sequence <- sample(bases, size=size, replace=TRUE)
  if (together) {
    sequence <- paste(sequence, collapse = "")
  }
  return(sequence)
}
```

We can fix this inability to generate multiple sequences by either
editing and adding to the functio body code (e.g. a for loop) or by
using the R **apply ** family of utility functions

``` r
sapply(6:12, generate_protein)
```

    [1] "AVATIA"       "NFIMSSR"      "CAIYFVAI"     "VEGSSWRKY"    "WQWEWQSRNH"  
    [6] "ETHIMNNREMG"  "GHNTWNVCVVHL"

``` r
generate_protein()
```

    [1] "ADIMYM"

> Q. Generate random protein sequences of length 6 to 12 amino acids.

``` r
ans <- sapply(6:12, generate_protein)

cat(ans, sep="\n")
```

    TDTKRR
    TKTSQRD
    YAMCTVMA
    CCIEDMCCQ
    QTHRQCGGMN
    NNTPDKVPKMV
    RCLERPPGDRFL

I want this to look like

    >ID.6
    EWIGVH
    >ID.7
    IHGPGHY
    >ID.8
    GQMKQNHS
    >ID.9
    WTVWYSTGI
    >ID.10
    TQSKKAGTRD
    >ID.11
    EVHTIMNENIG
    >ID.12
    FHYFDYYHIGAW

``` r
id_head <- function(number) {
  paste(">ID.", as.character(number), sep="")
}
proteins <- sapply(6:12, generate_protein)
heads <- sapply(6:12, id_head)
cat(paste(heads, proteins, sep="\n"), sep = "\n")
```

    >ID.6
    VNQGSV
    >ID.7
    QRMPIMC
    >ID.8
    ANWIRKAI
    >ID.9
    GCQTKTYGF
    >ID.10
    PDTSFLRKHW
    >ID.11
    NFGMHMGEDQT
    >ID.12
    FLIVDRLNGMGQ

``` r
cat(paste(">ID.", 6:12, "\n", ans, sep=""), sep="\n")
```

    >ID.6
    TDTKRR
    >ID.7
    TKTSQRD
    >ID.8
    YAMCTVMA
    >ID.9
    CCIEDMCCQ
    >ID.10
    QTHRQCGGMN
    >ID.11
    NNTPDKVPKMV
    >ID.12
    RCLERPPGDRFL

> Q. Determin if these sequences can be found in nature It would be cool
> and useful if I could get FASTA format output

I BLASTp searched my FASTA format sequences against NR and found that
length 6, 7 are not unique and can be found in the databases with 100%
coverage and 100% identity.

But for length 8, the coverage is 88% with 100% identity, or 100%
coverage and 87.50% identity. For length 9, 100% coverage and 89%
identity.

Random sequences of length 8 and above are unique and cannot be found in
database.

``` r
dim(cars)
```

    [1] 50  2
