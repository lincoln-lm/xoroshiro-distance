# xoroshiro-distance

C++ implementation of Pohlig-Hellman to solve the discrete logarithim in order to compute the "distance" (number of state transitions) between two states of xoroshiro128+. Factors 3, 5, 17, 257, and 641 are considered which allows accuracy in the range [0, 42007934]. Each smaller discrete logarithm is solved via lookup in a minimal perfect hash table.
