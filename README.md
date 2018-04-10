# hdsp
C++ implementation of DSP module for HTT allelic structure definition

Work in progress jesus christ i'm so rusty at C++

Utilises SEQAN for an alignment assembly interface (https://www.seqan.de/) so that must be in your /usr/local/include folder
Also utilises a C++ implementation of the difflib library from Python (https://github.com/duckie/difflib), pop that in /usr/local/include also.

Compiled via:
$ g++-7 -std=c++17 -lz -D SEQAN_HAS_ZLIB -fconcepts -I /usr/local/include
