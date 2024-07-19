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

discrete_distribution<int> build_frag_distribution(map<int, int> &frag_len_mapper, const char *frags_fp,
                                                   size_t seq_len, size_t ploidy,
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

    num_reads = std::round(((seq_len / avg_length) / ploidy) * coverage);

    return discrete_distribution<int>(frag_len_weights.begin(), frag_len_weights.end());
}

struct write_reads_args {
    string fasta_fp;
    string out_fp;
    string population;
    string ind;
    string chr;
    const map<int, int> *frag_len_mapper;
    discrete_distribution<int> *frag_dist;
    uniform_int_distribution<uint64_t> *pos_dist;
    size_t num_reads;
    size_t seq_len;
    mt19937 *generator;
};

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

void write_reads_with_errors(const write_reads_args &args)
{
    // Randomness engine
    std::random_device rd;
    mt19937 generator(rd());

    auto read_coords = simulate_read_coords(*(args.frag_len_mapper), *(args.frag_dist),
                                            *(args.pos_dist), args.num_reads,
                                            args.seq_len, generator);

    // TODO Test moving this to global
    uniform_real_distribution<float> error_dist(0.f, 1.f);
    uniform_int_distribution<int>    atgc_dist(0, 3);

    // Load true sequence
    auto true_sequence = load_sequence(args.fasta_fp, args.seq_len);

    // Write each read to FASTQ file
    std::ofstream out_f(args.out_fp, std::ofstream::out);
    string endo_read_tag = "@SEQ_" + args.population + "_" + args.ind + "_" + args.chr + "_";
    for (size_t i = 0; i != read_coords.size(); ++i) {
        uint64_t start = read_coords[i].first;
        uint64_t len   = read_coords[i].second;

        string read_tag = endo_read_tag + std::to_string(i);
        string quals(len, 'j');
        string to_write = read_tag + "\n";

        // Read that will get damaged
        string damaged = true_sequence.substr(start, start + len);

        vector<float> geno_errors(damaged.size());
        std::generate(geno_errors.begin(), geno_errors.end(),
                      [&](){return error_dist(generator);});

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

extern "C"
void simulate_reads(size_t fasta_fps_len, size_t seq_len, size_t ploidy, float coverage,
                    const char **fasta_fps, const char **populations,
                    const char **individuals, const char **chromosomes, const char **out_fps,
                    const char *cont_fp, const char *frags_fp)
{
    // Build generator for fragment lengths and starts
    map<int, int> frag_len_mapper;
    size_t num_reads = 0;
    auto frag_dist = build_frag_distribution(frag_len_mapper, frags_fp, seq_len,
                                             ploidy, coverage, num_reads);
    uniform_int_distribution<uint64_t> pos_dist(0, seq_len - 1);


    // TODO Handle contamination sequence
    // TODO Handle misincorporation stuff

    #pragma omp parallel for
    for (size_t i = 0; i != fasta_fps_len; ++i) {
        write_reads_args new_args;
        new_args.frag_len_mapper = &frag_len_mapper;
        new_args.fasta_fp   = fasta_fps[i];
        new_args.out_fp     = out_fps[i];
        new_args.population = populations[i];
        new_args.ind        = individuals[i];
        new_args.chr        = chromosomes[i];
        new_args.frag_dist  = &frag_dist;
        new_args.pos_dist   = &pos_dist;
        new_args.num_reads  = num_reads;
        new_args.seq_len    = seq_len;
        //new_args.generator  = &generator;
        write_reads_with_errors(new_args);
    }

    return;
}
