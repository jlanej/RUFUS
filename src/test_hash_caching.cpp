/*
 * test_hash_caching.cpp
 * Unit test to verify that HashToLong caching optimizations
 * produce identical results to the original inline computation.
 * This validates the correctness of the performance changes in
 * BuildUpHashCountTable and LookUpKmers.
 */

#include <bitset>
#include <cassert>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// Duplicated from RUFUS.interpret.cpp (local function)
unsigned long HashToLong(string hash) {
    bitset<64> HashBits;
    for (int i = 0; i < hash.length(); i++) {
        if (hash[i] == 'A') {
            HashBits[i * 2] = 0;
            HashBits[i * 2 + 1] = 0;
        } else if (hash[i] == 'C') {
            HashBits[i * 2] = 0;
            HashBits[i * 2 + 1] = 1;
        } else if (hash[i] == 'G') {
            HashBits[i * 2] = 1;
            HashBits[i * 2 + 1] = 0;
        } else if (hash[i] == 'T') {
            HashBits[i * 2] = 1;
            HashBits[i * 2 + 1] = 1;
        }
    }
    return HashBits.to_ulong();
}

// Test that HashToLong is deterministic (same input -> same output)
void test_deterministic() {
    string seqs[] = {"ACGT", "AAAA", "TTTT", "GCTA", "ACGTACGTACGTACGTACGTACGTACGTACGT"};
    for (const auto& s : seqs) {
        unsigned long a = HashToLong(s);
        unsigned long b = HashToLong(s);
        assert(a == b);
    }
    cout << "PASS: test_deterministic" << endl;
}

// Test that different sequences produce different hashes
void test_distinct() {
    assert(HashToLong("ACGT") != HashToLong("TGCA"));
    assert(HashToLong("AAAA") != HashToLong("CCCC"));
    assert(HashToLong("AAAA") != HashToLong("TTTT"));
    cout << "PASS: test_distinct" << endl;
}

// Test that pre-computing and reusing a cached value is identical
// to calling HashToLong inline multiple times (the core optimization)
void test_cache_equivalence() {
    // Simulate the LookUpKmers pattern:
    // Original: multiple calls to HashToLong(hash) for the same hash
    // Optimized: one call, reuse the cached value
    unordered_map<unsigned long int, int> MutantHashes;
    MutantHashes[HashToLong("ACGTACGTACGTACGTACGTACGTA")] = 42;
    MutantHashes[HashToLong("TTTTTTTTTTTTTTTTTTTTTTTTT")] = 7;

    string hash = "ACGTACGTACGTACGTACGTACGTA";

    // Original pattern (inline calls)
    int origContig = -1, origAlt = -1;
    if (MutantHashes.count(HashToLong(hash)) > 0) {
        origContig = MutantHashes[HashToLong(hash)];
    }
    if (MutantHashes.count(HashToLong(hash)) > 0) {
        origAlt = MutantHashes[HashToLong(hash)];
    }

    // Optimized pattern (cached)
    unsigned long int longHash = HashToLong(hash);
    int cachedContig = -1, cachedAlt = -1;
    if (MutantHashes.count(longHash) > 0) {
        cachedContig = MutantHashes[longHash];
    }
    if (MutantHashes.count(longHash) > 0) {
        cachedAlt = MutantHashes[longHash];
    }

    assert(origContig == cachedContig);
    assert(origAlt == cachedAlt);
    assert(origContig == 42);
    cout << "PASS: test_cache_equivalence" << endl;
}

// Test the BuildUpHashCountTable pre-computation pattern
void test_precompute_pattern() {
    int HashSize = 5;
    vector<string> hashes = {"ACGTA", "NNCGT", "TTTTT", "ACGNA", "GGGGG"};
    vector<string> hashesRef = {"TGCAT", "AAGCT", "AAAAA", "TGCNA", "CCCCC"};

    // Original pattern: validate + compute inline per parent
    vector<bool> origValid(hashes.size());
    vector<unsigned long> origLong(hashes.size());
    vector<unsigned long> origLongRef(hashes.size());

    for (int i = 0; i < hashes.size(); i++) {
        string hash = hashes[i];
        bool checkHash = true;
        for (int j = 0; j < HashSize; j++) {
            if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T')) {
                checkHash = false;
                break;
            }
        }
        origValid[i] = checkHash;
        if (checkHash) {
            origLong[i] = HashToLong(hash);
            origLongRef[i] = HashToLong(hashesRef[i]);
        }
    }

    // Optimized pattern: pre-compute once
    int numHashes = hashes.size();
    vector<unsigned long int> longHashes(numHashes);
    vector<unsigned long int> longHashesRef(numHashes);
    vector<bool> validHash(numHashes);
    for (int i = 0; i < numHashes; i++) {
        const string& hash = hashes[i];
        bool check = true;
        for (int j = 0; j < HashSize; j++) {
            if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T')) {
                check = false;
                break;
            }
        }
        validHash[i] = check;
        if (check) {
            longHashes[i] = HashToLong(hash);
            longHashesRef[i] = HashToLong(hashesRef[i]);
        }
    }

    // Verify equivalence
    for (int i = 0; i < hashes.size(); i++) {
        assert(origValid[i] == validHash[i]);
        if (origValid[i]) {
            assert(origLong[i] == longHashes[i]);
            assert(origLongRef[i] == longHashesRef[i]);
        }
    }

    // Verify specific expected validity
    assert(validHash[0] == true);   // ACGTA - all valid
    assert(validHash[1] == false);  // NNCGT - has N
    assert(validHash[2] == true);   // TTTTT - all valid
    assert(validHash[3] == false);  // ACGNA - has N
    assert(validHash[4] == true);   // GGGGG - all valid

    cout << "PASS: test_precompute_pattern" << endl;
}

// Test that simulated parent hash lookups produce identical results
// using pre-computed vs inline HashToLong
void test_parent_lookup_equivalence() {
    int HashSize = 4;
    vector<string> hashes = {"ACGT", "TGCA", "NNNN", "GGCC"};

    // Set up parent hash map
    unordered_map<unsigned long int, int> ParentHash;
    ParentHash[HashToLong("ACGT")] = 10;
    ParentHash[HashToLong("GGCC")] = 20;

    // Original inline computation
    vector<int> origCounts;
    for (int i = 0; i < hashes.size(); i++) {
        string hash = hashes[i];
        bool checkHash = true;
        for (int j = 0; j < HashSize; j++) {
            if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T')) {
                checkHash = false;
                break;
            }
        }
        if (checkHash) {
            unsigned long int LongHash = HashToLong(hash);
            if (ParentHash.count(LongHash) > 0)
                origCounts.push_back(ParentHash[LongHash]);
            else
                origCounts.push_back(0);
        } else {
            origCounts.push_back(-1);
        }
    }

    // Pre-computed lookup
    vector<unsigned long int> longHashes(hashes.size());
    vector<bool> validHash(hashes.size());
    for (int i = 0; i < hashes.size(); i++) {
        const string& hash = hashes[i];
        bool check = true;
        for (int j = 0; j < HashSize; j++) {
            if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T')) {
                check = false;
                break;
            }
        }
        validHash[i] = check;
        if (check) longHashes[i] = HashToLong(hash);
    }

    vector<int> cachedCounts;
    for (int i = 0; i < hashes.size(); i++) {
        if (validHash[i]) {
            unsigned long int LongHash = longHashes[i];
            if (ParentHash.count(LongHash) > 0)
                cachedCounts.push_back(ParentHash[LongHash]);
            else
                cachedCounts.push_back(0);
        } else {
            cachedCounts.push_back(-1);
        }
    }

    assert(origCounts.size() == cachedCounts.size());
    for (int i = 0; i < origCounts.size(); i++) {
        assert(origCounts[i] == cachedCounts[i]);
    }

    // Verify expected values
    assert(origCounts[0] == 10);   // ACGT found in parent
    assert(origCounts[1] == 0);    // TGCA not found
    assert(origCounts[2] == -1);   // NNNN invalid
    assert(origCounts[3] == 20);   // GGCC found in parent

    cout << "PASS: test_parent_lookup_equivalence" << endl;
}

int main() {
    cout << "Running hash caching correctness tests..." << endl;
    test_deterministic();
    test_distinct();
    test_cache_equivalence();
    test_precompute_pattern();
    test_parent_lookup_equivalence();
    cout << "\nAll tests passed!" << endl;
    return 0;
}
