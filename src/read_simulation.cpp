#include <iostream>
#include <random>
#include <vector>
#include <map>
#include <cmath>

using std::vector;
using std::map;

std::discrete_distribution<int> build_frag_distribution(map<int, uint64_t> &frag_len_mapper, const char *frags_fp,
                                                        uint64_t seq_len, size_t ploidy,
                                                        double coverage, size_t &num_reads)
{
    float avg_length = 0;
    vector<double> frag_len_weights;
    if (frags_fp[0] != '\0') {
        // TODO Handle frag distribution file
    } else {
        frag_len_mapper[0] = 70;
        frag_len_weights.push_back(1);
        avg_length = 70;
    }

    num_reads = std::round(((seq_len / avg_length) / ploidy) * coverage)

    return std::discrete_distribution<int>(frag_len_weights.begin(), frag_len_weights.end());
}

extern "C"
void simulate_reads(size_t fasta_fps_len, const char **fasta_fps, const char *cont_fp,
                    const char *frags_fp, size_t seq_len, size_t ploidy, double coverage)
{
    // Build generator for fragment lengths
    map<int, uint64_t> frag_len_mapper;
    size_t num_reads = 0;
    auto frag_distribution = build_frag_distribution(frag_len_mapper, frags_fp, seq_len,
                                                     ploidy, coverage, num_reads);

    // TODO Handle contamination sequence
    // TODO Handle misincorporation stuff
    return;
}
