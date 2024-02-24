# Algorithms-MinHash
This program use the Min Hash algorithm to compute the Jaccard Similarity between files.

### Implementation
1. The files are shingled first to save space.
2. Using upper triangular matrix instead of plain matrix to save space.
3. Implement write buffer, small writes are firstly stored in the memory, only when the buffer becomes full that it writes to the disk. This reduce the number of write, thus save running time.
4. Using Locality-Sensitive Hasing to further simplify the comparison between signature matrix.

### Files
1. min_hash.c min_hash.h		the min_hash algorithm, generating files automatically and print the result.
2. CRC32.cpp CRC32.h CRC32 hash algorithm that is used in min hash
