#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>
#include <execution>
#include <utility>

using std::vector;      using std::pair;
using std::map;         using std::mt19937;
using std::string;
using std::discrete_distribution;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

////////////////
// STRUCTURES //
////////////////

// Simple structure to store misincorporation profile data
struct misincorporation_pkg {
    size_t max_pos;
};

///////////////
// FUNCTIONS //
///////////////

// Parse fragmentation distribution file into a discrete distribution
discrete_distribution<int> build_frag_distribution(map<int, int> &frag_len_mapper, size_t frag_dist_len,
                                                   const int *frag_lens, const double *frag_probs)
{
    for (size_t i = 0; i != frag_dist_len; ++i)
        frag_len_mapper[i] = frag_lens[i];

    vector<double> frag_dist_weights(frag_probs, frag_probs + frag_dist_len);
    return discrete_distribution<int>(frag_dist_weights.begin(), frag_dist_weights.end());
}

// Load FASTA file into memory
string load_sequence(const string &fp, size_t seq_len)
{
    string         sequence;
    string         line;
    std::ifstream  f(fp);

    sequence.reserve(seq_len);
    while (std::getline(f, line)) {
        if (line[0] == '>')
            continue;
        sequence.append(line);
    }

    return sequence;
}

// Simulate starts and lengths of read fragments
vector<pair<uint64_t, uint64_t>> simulate_read_coords(const map<int, int> &frag_len_mapper,
                                                      discrete_distribution<int> &frag_dist,
                                                      uniform_int_distribution<uint64_t> &pos_dist,
                                                      size_t num_reads, size_t seq_len, mt19937 &generator)
{
    vector<pair<uint64_t, uint64_t>> ret(num_reads);

    for (auto &p : ret) {
        uint64_t start = pos_dist(generator);
        uint64_t len   = frag_len_mapper.at(frag_dist(generator));
        p.first  = start;
        p.second = start + len - 1 >= seq_len ? seq_len - start : len;
    }

    return ret;
}

// Main function of interest, simulate reads with genotyping error and misincorporations
void write_reads_with_errors(const string &fasta_fp, const string &out_fp, const string &population,
                             const string &ind, const string &chr, const map<int, int> &frag_len_mapper,
                             discrete_distribution<int> &frag_dist, size_t num_reads, size_t seq_len,
                             double cont_p, bool has_cont, const string &cont_seq)
{
    // Randomness engine
    std::random_device rd;
    mt19937 generator(rd());

    // Distributions for randomness. TODO Move to shared memory?
    uniform_real_distribution<double>  error_dist(0.f, 1.f);
    uniform_int_distribution<int>      atgc_dist(0, 3);
    uniform_int_distribution<uint64_t> pos_dist(0, seq_len - 1);

    // Simulate positions of reads
    auto read_coords = simulate_read_coords(frag_len_mapper, frag_dist, pos_dist, num_reads, seq_len, generator);

    // Load true sequence
    auto true_sequence = load_sequence(fasta_fp, seq_len);

    // Write each read to FASTQ file
    std::ofstream out_f(out_fp, std::ofstream::out);
    string endo_read_tag = "@SEQ_" + population + "_" + ind + "_" + chr + "_";
    string cont_read_tag = "@SEQ_CONTAMINATED_" + population + "_" + ind + "_" + chr + "_";
    for (size_t i = 0; i != read_coords.size(); ++i) {
        uint64_t start = read_coords[i].first;
        uint64_t len   = read_coords[i].second;
        string quals(len, 'j');

        // Introduce contamination if contamination sequence was provided
        if (has_cont && error_dist(generator) < cont_p) {
            out_f << cont_read_tag + std::to_string(i) + "\n" + cont_seq.substr(start, len) + "\n" + quals + "\n";
            continue;
        }

        string read_tag = endo_read_tag + std::to_string(i);
        string to_write = read_tag + "\n";

        // Read that will get damaged
        string damaged = true_sequence.substr(start, len);

        // Generate dice rolls for genotyping error at each base
        vector<double> geno_errors(damaged.size());
        std::generate(geno_errors.begin(), geno_errors.end(),
                      [&](){return error_dist(generator);});

        // Check each base for genotyping and misincorporation damage
        for (size_t j = 0; j != damaged.size(); ++j) {
            if (geno_errors[j] < 0.001333333) {
                damaged[j] = "ATGC"[atgc_dist(generator)];
            }
        }

        // Write
        to_write += damaged + "\n" + quals + "\n";
        out_f << to_write;
    }
}

// Entry point for Python data
extern "C"
void simulate_reads(size_t fasta_fps_len, size_t seq_len, size_t num_reads, size_t frag_dist_len,
                    const int *frag_lens, const double *frag_probs, double cont_p,
                    const char **fasta_fps, const char **populations, const char **individuals,
                    const char **chromosomes, const char **out_fps, const char *cont_fp)
{
    // Build generator for fragment lengths and starts
    map<int, int> frag_len_mapper;
    auto frag_dist = build_frag_distribution(frag_len_mapper, frag_dist_len, frag_lens, frag_probs);

    // Handle contamination sequence
    string cont_sequence = "";
    bool   has_cont = (cont_fp[0] != '\0');
    if (has_cont)
        cont_sequence = load_sequence(cont_fp, seq_len);

    // TODO Handle misincorporation stuff

    #pragma omp parallel for
    for (size_t i = 0; i != fasta_fps_len; ++i) {
        write_reads_with_errors(fasta_fps[i], out_fps[i], populations[i], individuals[i], chromosomes[i],
                                frag_len_mapper, frag_dist, num_reads, seq_len, cont_p, has_cont, cont_sequence);
    }
}
