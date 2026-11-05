// Exhaustive search for a kernel solution that uses exactly K centres,
// with replacement, supported on a given eligible column list.
#include "kernel_common.hpp"
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "usage: enum_k N r m K eligible.txt\n";
        return 2;
    }
    int N = std::atoi(argv[1]), r = std::atoi(argv[2]), m = std::atoi(argv[3]);
    int K = std::atoi(argv[4]);
    std::ifstream in(argv[5]);
    std::vector<int> elig;
    int x;
    while (in >> x) elig.push_back(x);
    auto T = build_kernel(N, r, m);
    const int nrows = T.n_pts;
    const int E = static_cast<int>(elig.size());
    std::cerr << "K " << K << " eligible " << E << " points " << nrows << "\n";
    std::vector<int64_t> cols(static_cast<size_t>(E) * nrows);
    for (int i = 0; i < E; ++i)
        std::copy(T.scaled.data() + static_cast<size_t>(elig[i]) * nrows,
                  T.scaled.data() + static_cast<size_t>(elig[i] + 1) * nrows,
                  cols.begin() + static_cast<size_t>(i) * nrows);
    std::vector<int64_t> row_max(nrows, 0);
    for (int i = 0; i < E; ++i) {
        const int64_t* c = cols.data() + static_cast<size_t>(i) * nrows;
        for (int row = 0; row < nrows; ++row) row_max[row] = std::max(row_max[row], c[row]);
    }
    // Any partial sum that cannot reach the rhs even with the remaining
    // centres at row_max is dead. Order rows by a short random sample.
    std::vector<int> fail(nrows, 0);
    std::mt19937 rng(1);
    std::uniform_int_distribution<int> dist(0, E - 1);
    std::vector<int> pick(K);
    for (int s = 0; s < 2000; ++s) {
        for (int t = 0; t < K; ++t) pick[t] = dist(rng);
        for (int row = 0; row < nrows; ++row) {
            int64_t ssum = 0;
            for (int t = 0; t < K; ++t)
                ssum += cols[static_cast<size_t>(pick[t]) * nrows + row];
            if (ssum < T.rhs[row]) { fail[row]++; break; }
        }
    }
    std::vector<int> order(nrows);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return fail[a] > fail[b]; });

    std::atomic<int64_t> solutions{0};
    std::atomic<int> next_a{0};
    auto t0 = std::chrono::steady_clock::now();
    unsigned nthreads = 4u;

    auto worker = [&]() {
        std::vector<int64_t> acc(static_cast<size_t>(K + 1) * nrows, 0);
        std::vector<int> choice(K);
        auto dfs = [&](auto&& self, int depth, int start) -> void {
            if (solutions.load()) return;
            if (depth == K) {
                const int64_t* s = acc.data() + static_cast<size_t>(depth) * nrows;
                for (int row : order) if (s[row] < T.rhs[row]) return;
                solutions.fetch_add(1);
                std::cerr << "SOL";
                for (int t = 0; t < K; ++t) std::cerr << " " << elig[choice[t]];
                std::cerr << "\n";
                return;
            }
            const int64_t* prev = acc.data() + static_cast<size_t>(depth) * nrows;
            int64_t* cur = acc.data() + static_cast<size_t>(depth + 1) * nrows;
            int left = K - depth - 1;
            for (int i = start; i < E; ++i) {
                if (solutions.load()) return;
                const int64_t* col = cols.data() + static_cast<size_t>(i) * nrows;
                bool ok = true;
                for (int row : order) {
                    __int128 v = static_cast<__int128>(prev[row]) + col[row]
                                 + static_cast<__int128>(left) * row_max[row];
                    if (v < T.rhs[row]) { ok = false; break; }
                    cur[row] = prev[row] + col[row];
                }
                if (!ok) continue;
                // fill any rows not written because we broke? we write only until break.
                // Recompute cur fully. The break may have left cur partial.
                for (int row = 0; row < nrows; ++row) cur[row] = prev[row] + col[row];
                choice[depth] = i;
                self(self, depth + 1, i);
            }
        };
        while (true) {
            int a = next_a.fetch_add(1);
            if (a >= E) return;
            const int64_t* col = cols.data() + static_cast<size_t>(a) * nrows;
            int64_t* cur = acc.data() + nrows;
            for (int row = 0; row < nrows; ++row) cur[row] = col[row];
            choice[0] = a;
            dfs(dfs, 1, a);
            if (a % 5 == 0) {
                auto sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                std::cerr << "a " << a << " sol " << solutions.load() << " sec " << sec << "\n";
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nthreads; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    auto sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "solutions " << solutions.load() << " seconds " << sec << "\n";
    return solutions.load() ? 1 : 0;
}
