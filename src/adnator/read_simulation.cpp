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
    bool valid;
    map<pair<size_t, char>, discrete_distribution<int>> mis_prob;
};

///////////////
// FUNCTIONS //
///////////////

// Turn ctypes friendly version of the data into a misincorporation package struct
misincorporation_pkg build_mis_mapper(size_t mis_len, const size_t *mis_pos,
                                      const char *mis_nuc, const double *mis_pro)
{
    misincorporation_pkg ret;

    if (mis_len == 0) {
        ret.valid = false;
        return ret;
    }

    ret.valid = true;
    ret.max_pos = 0;

    for (size_t i = 0; i != mis_len; ++i) {
        ret.max_pos = std::max(ret.max_pos, mis_pos[i]);
        vector<double> probs(mis_pro + (i * 4), mis_pro + (i * 4) + 4);
        ret.mis_prob[std::make_pair(mis_pos[i], mis_nuc[i])] =
            discrete_distribution<int>(probs.begin(), probs.end());
    }

    return ret;
}

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
                             double cont_p, bool has_cont, const string &cont_seq,
                             misincorporation_pkg &mis_5_mapper, misincorporation_pkg &mis_3_mapper,
                             bool has_geno_err, int random_seed)
{
    // Randomness engine
    mt19937 generator(random_seed);

    // Distributions for randomness. TODO Move to shared memory?
    uniform_real_distribution<double>  error_dist(0.f, 1.f);
    uniform_int_distribution<int>      atgc_dist(0, 3);
    uniform_int_distribution<uint64_t> pos_dist(0, seq_len - 1);

    // Simulate positions of reads
    auto read_coords = simulate_read_coords(frag_len_mapper, frag_dist, pos_dist, num_reads, seq_len, generator);

    // Load true sequence
    auto true_sequence = load_sequence(fasta_fp, seq_len);

    // Write each read to FASTQ file
    std::ofstream out_f(out_fp, std::ios::app);
    string read_tag_info = population + "_" + ind + "_" + chr + "_";
    string endo_read_tag = "@SEQ_" + read_tag_info;
    string cont_read_tag = "@SEQ_CONTAMINATED_" + read_tag_info;
    for (size_t i = 0; i != read_coords.size(); ++i) {
        uint64_t start = read_coords[i].first;
        uint64_t len   = read_coords[i].second;
        string quals(len, 'j');
        string curr_read_info = std::to_string(i) + "_" + std::to_string(start);

        // Introduce contamination if contamination sequence was provided
        if (has_cont && error_dist(generator) < cont_p) {
            out_f << cont_read_tag + curr_read_info + "\n" + cont_seq.substr(start, len) + "\n" + quals + "\n";
            continue;
        }

        string to_write = endo_read_tag + curr_read_info + "\n";

        // Read that will get damaged
        string damaged = true_sequence.substr(start, len);

        // Generate dice rolls for genotyping error at each base
        vector<double> geno_errors(damaged.size());
        if (has_geno_err) {
            std::generate(geno_errors.begin(), geno_errors.end(),
                          [&](){return error_dist(generator);});
        }

        // Check each base for misincorporation damage and genotyping error
        for (size_t j = 0; j != damaged.size(); ++j) {
            if (damaged[j] != 'A' && damaged[j] != 'T' && damaged[j] != 'G' && damaged[j] != 'C')
                continue;

            if (mis_5_mapper.valid && j <= mis_5_mapper.max_pos) {
                damaged[j] = "ATGC"[mis_5_mapper.mis_prob.at({j, damaged[j]})(generator)];
            } else if (mis_3_mapper.valid && (len - (j + 1)) <= mis_3_mapper.max_pos) {
                damaged[j] = "ATGC"[mis_3_mapper.mis_prob.at({(len - (j + 1)), damaged[j]})(generator)];
            }

            if (has_geno_err && geno_errors[j] < 0.001333333) {
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
                    const char **chromosomes, const char **out_fps, const char *cont_fp,
                    size_t mis_5_len, const size_t *mis_5_pos, const char *mis_5_nuc, const double *mis_5_pro,
                    int mis_3_len, const size_t *mis_3_pos, const char *mis_3_nuc, const double *mis_3_pro,
                    bool has_geno_err, int random_seed, size_t ploidy)
{
    // Build generator for fragment lengths and starts
    map<int, int> frag_len_mapper;
    auto frag_dist = build_frag_distribution(frag_len_mapper, frag_dist_len, frag_lens, frag_probs);

    // Build mappers of position and nucleotide to misincorporation probabilities
    auto mis_5_mapper = build_mis_mapper(mis_5_len, mis_5_pos, mis_5_nuc, mis_5_pro);
    auto mis_3_mapper = build_mis_mapper(mis_3_len, mis_3_pos, mis_3_nuc, mis_3_pro);

    // Handle contamination sequence
    string cont_sequence = "";
    bool   has_cont = (cont_fp[0] != '\0');
    if (has_cont)
        cont_sequence = load_sequence(cont_fp, seq_len);

    // Generate random seeds for all parallelized iterations
    mt19937 main_rng;
    vector<int> random_seeds;
    if (random_seed == -1) {
        std::random_device rd;
        main_rng.seed(rd());
    } else {
        main_rng.seed(random_seed);
    }

    for (size_t i = 0; i != fasta_fps_len; ++i)
        random_seeds.push_back(main_rng());

    # pragma omp parallel for
    for (size_t i = 0; i < fasta_fps_len; i += ploidy) {
        for (size_t j = i; j != i + ploidy; ++j) {
            write_reads_with_errors(fasta_fps[j], out_fps[j], populations[j], individuals[j], chromosomes[j],
                                    frag_len_mapper, frag_dist, num_reads, seq_len, cont_p, has_cont, cont_sequence,
                                    mis_5_mapper, mis_3_mapper, has_geno_err, random_seeds[j]);
        }
    }
}
