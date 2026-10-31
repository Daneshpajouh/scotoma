// Exhaustive search for nonnegative integer class-multiplicity solutions
// of sum <= KMAX. Prints every solution with sum == KMAX (and smaller ones).
// Exact int64 arithmetic. Parallel over the first index, one thread.
#include "kernel_common.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <mutex>
#include <numeric>
#include <random>
#include <thread>
#include <vector>

static bool covers_sum3(const int64_t* a, const int64_t* b, const int64_t* c,
                        const int64_t* rhs, const int* order, int nrows) {
    for (int t = 0; t < nrows; ++t) {
        int r = order[t];
        if (a[r] + b[r] + c[r] < rhs[r]) return false;
    }
    return true;
}

static bool covers_sum2(const int64_t* a, const int64_t* b,
                        const int64_t* rhs, const int* order, int nrows) {
    for (int t = 0; t < nrows; ++t) {
        int r = order[t];
        if (a[r] + b[r] < rhs[r]) return false;
    }
    return true;
}

static bool covers_sum1(const int64_t* a, const int64_t* rhs, const int* order, int nrows) {
    for (int t = 0; t < nrows; ++t) {
        int r = order[t];
        if (a[r] < rhs[r]) return false;
    }
    return true;
}

int main(int argc, char** argv) {
    int N = 19, r = 4, m = 13, kmax = 3;
    std::string out_path = "certificates/class_triples_19_4_13.txt";
    if (argc >= 4) {
        N = std::atoi(argv[1]);
        r = std::atoi(argv[2]);
        m = std::atoi(argv[3]);
    }
    if (argc >= 5) kmax = std::atoi(argv[4]);
    if (argc >= 6) out_path = argv[5];
    if (kmax != 3) {
        std::cerr << "this enumerator handles kmax=3 only\n";
        return 2;
    }
    auto T = build_kernel(N, r, m);
    const int n = T.n_cent;
    const int nrows = T.n_pts;
    std::cerr << "cell " << N << " " << r << " " << m
              << " centres " << n << " points " << nrows << "\n";

    // Order rows by a sample: rows that reject random triples earliest go first.
    std::vector<int> fail_count(nrows, 0);
    std::mt19937 rng(1);
    std::uniform_int_distribution<int> dist(0, n - 1);
    // mt19937 needs include
    for (int s = 0; s < 2000; ++s) {
        int i = dist(rng), j = dist(rng), k = dist(rng);
        const int64_t* a = T.scaled.data() + static_cast<size_t>(i) * nrows;
        const int64_t* b = T.scaled.data() + static_cast<size_t>(j) * nrows;
        const int64_t* c = T.scaled.data() + static_cast<size_t>(k) * nrows;
        for (int row = 0; row < nrows; ++row) {
            if (a[row] + b[row] + c[row] < T.rhs[row]) {
                fail_count[row]++;
                break;
            }
        }
    }
    std::vector<int> order(nrows);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return fail_count[a] > fail_count[b];
    });

    std::vector<int64_t> row_max(nrows, 0);
    for (int c = 0; c < n; ++c) {
        const int64_t* col = T.scaled.data() + static_cast<size_t>(c) * nrows;
        for (int row = 0; row < nrows; ++row) row_max[row] = std::max(row_max[row], col[row]);
    }

    std::atomic<int64_t> n1{0}, n2{0}, n3{0}, seen{0}, pairs_kept{0};
    std::mutex file_mu;
    std::ofstream out(out_path);
    if (!out) {
        std::cerr << "cannot write " << out_path << "\n";
        return 1;
    }
    out << "N " << N << " r " << r << " m " << m << " kmax " << kmax
        << " centres " << n << " points " << nrows << "\n";
    out << "format: sum i j k then prefix:weight of each class, i<=j<=k\n";
    // Keep every solution up to a cap, then count only.
    const int64_t write_cap = 20000;

    auto dump = [&](int i, int j, int k) {
        if (n3.load() >= write_cap) return;
        std::lock_guard<std::mutex> lock(file_mu);
        if (n3.load() >= write_cap) return;
        auto emit = [&](int c) {
            int pref = c / (T.shape.C + 1);
            int w = c % (T.shape.C + 1);
            out << " " << prefix_string(pref, T.shape.L) << ":" << w;
        };
        out << "3 " << i << " " << j << " " << k;
        emit(i); emit(j); emit(k);
        out << "\n";
    };

    unsigned nthreads = 1u;
    if (nthreads < 1) nthreads = 1;
    std::cerr << "threads " << nthreads << "\n";
    std::atomic<int> next_i{0};
    auto worker = [&]() {
        while (true) {
            int i = next_i.fetch_add(1);
            if (i >= n) return;
            const int64_t* ai = T.scaled.data() + static_cast<size_t>(i) * nrows;
            if (covers_sum1(ai, T.rhs.data(), order.data(), nrows)) {
                // 3*column: check separately below; 1*column is a 1-cover
                n1.fetch_add(1);
            }
            // triple with the same column three times, and pairs
            std::vector<int64_t> triple(nrows), doub(nrows);
            for (int row = 0; row < nrows; ++row) {
                doub[row] = ai[row] + ai[row];
                triple[row] = doub[row] + ai[row];
            }
            if (covers_sum1(triple.data(), T.rhs.data(), order.data(), nrows)) {
                // recorded as i,i,i when we reach j=k=i
            }
            for (int j = i; j < n; ++j) {
                const int64_t* aj = T.scaled.data() + static_cast<size_t>(j) * nrows;
                bool pair_ok = true;
                bool extendable = true;
                for (int t = 0; t < nrows; ++t) {
                    int row = order[t];
                    int64_t s2 = ai[row] + aj[row];
                    if (s2 < T.rhs[row]) pair_ok = false;
                    if (s2 + row_max[row] < T.rhs[row]) { extendable = false; break; }
                }
                if (!extendable) continue;
                pairs_kept.fetch_add(1);
                if (pair_ok) n2.fetch_add(1);
                for (int k = j; k < n; ++k) {
                    const int64_t* ak = T.scaled.data() + static_cast<size_t>(k) * nrows;
                    if (covers_sum3(ai, aj, ak, T.rhs.data(), order.data(), nrows)) {
                        dump(i, j, k);
                        n3.fetch_add(1);
                    }
                }
            }
            int64_t s = seen.fetch_add(1) + 1;
            if (s % 25 == 0) {
                std::cerr << "i_done " << s << " solutions " << n3.load() << "\n";
            }
        }
    };
    auto t0 = std::chrono::steady_clock::now();
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nthreads; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    auto sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out << "counts 1covers " << n1.load() << " 2covers " << n2.load()
        << " 3covers " << n3.load() << " seconds " << sec << "\n";
    out.flush();
    std::cerr << "done 1=" << n1.load() << " 2=" << n2.load() << " 3=" << n3.load()
              << " pairs_kept " << pairs_kept.load() << " sec " << sec << "\n";
    std::cout << "SOLUTIONS " << n3.load() << "\n";
    return 0;
}
