// Exact class-kernel arithmetic for the canonical single-deletion cell.
// A centre class is a prefix in {0,1,2}^L together with a suffix weight.
// Letter 2 stands for either of the two letters the library never uses.
// Point classes are the same kind of pair, restricted to those inside U.
//
// scaled(c, p) / rhs(p) equals the fraction of point-class p covered by any
// one concrete centre of class c. Comparison of integer sums against rhs is
// exact. No floating point is used.
#pragma once
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct CellShape {
    int N = 0, r = 0, m = 0, L = 0, C = 0;
};

inline CellShape make_shape(int N, int r, int m) {
    if (m < 0 || m > N || r < 1 || r > N + 1) throw std::runtime_error("bad cell");
    CellShape s;
    s.N = N;
    s.r = r;
    s.m = m;
    s.L = r - 1;
    s.C = N - r + 1;
    if (s.L < 0 || s.C < 0) throw std::runtime_error("bad shape");
    return s;
}

inline int binom_small(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    int64_t v = 1;
    for (int i = 1; i <= k; ++i) v = v * (n - k + i) / i;
    return static_cast<int>(v);
}

inline int ipow_int(int base, int exp) {
    int v = 1;
    for (int i = 0; i < exp; ++i) v *= base;
    return v;
}

// prefix id in lexicographic order, digits in {0,1,2}, length L
inline void id_to_prefix(int id, int L, int* out) {
    for (int i = L - 1; i >= 0; --i) {
        out[i] = id % 3;
        id /= 3;
    }
}

inline int prefix_to_id(const int* a, int L) {
    int id = 0;
    for (int i = 0; i < L; ++i) id = id * 3 + a[i];
    return id;
}

struct KernelTable {
    CellShape shape;
    int n_pref = 0;
    int n_cent = 0;
    int n_pts = 0;
    std::vector<int> point_pref;   // size n_pts
    std::vector<int> point_w;      // size n_pts
    std::vector<int64_t> rhs;      // size n_pts
    // column-major scaled numerators, index centre * n_pts + point
    std::vector<int64_t> scaled;
    std::vector<uint8_t> dmin;     // per prefix id

    int64_t at(int centre, int point) const { return scaled[static_cast<size_t>(centre) * n_pts + point]; }
};

inline int64_t pow3_table[32];
inline int64_t pow2_table[32];

inline void init_pows() {
    pow3_table[0] = 1;
    pow2_table[0] = 1;
    for (int i = 1; i < 32; ++i) {
        pow3_table[i] = pow3_table[i - 1] * 3;
        pow2_table[i] = pow2_table[i - 1] * 2;
    }
}

// counts[p * (C+1) * (C+1) + w * (C+1) + d]
inline std::vector<int64_t> suffix_counts(int C) {
    const int W = C + 1;
    std::vector<int64_t> counts(static_cast<size_t>(W) * W * W, 0);
    std::vector<int64_t> cur(static_cast<size_t>(W) * W, 0), nxt(static_cast<size_t>(W) * W, 0);
    for (int p = 0; p <= C; ++p) {
        std::fill(cur.begin(), cur.end(), 0);
        cur[0] = 1;
        for (int coord = 0; coord < C; ++coord) {
            std::fill(nxt.begin(), nxt.end(), 0);
            const bool hot = coord < p;
            for (int w = 0; w <= C; ++w) {
                for (int d = 0; d <= C; ++d) {
                    int64_t n = cur[w * W + d];
                    if (!n) continue;
                    if (hot) {
                        // point 0: weight 0, distance 1, one letter
                        if (d + 1 <= C) nxt[w * W + (d + 1)] += n;
                        // point equals centre: weight 1, distance 0, one letter
                        if (w + 1 <= C) nxt[(w + 1) * W + d] += n;
                        // other nonzero: weight 1, distance 1, two letters
                        if (w + 1 <= C && d + 1 <= C) nxt[(w + 1) * W + (d + 1)] += 2 * n;
                    } else {
                        // point 0: weight 0, distance 0
                        nxt[w * W + d] += n;
                        // nonzero: weight 1, distance 1, three letters
                        if (w + 1 <= C && d + 1 <= C) nxt[(w + 1) * W + (d + 1)] += 3 * n;
                    }
                }
            }
            cur.swap(nxt);
        }
        for (int w = 0; w <= C; ++w) {
            int64_t sum = 0;
            for (int d = 0; d <= C; ++d) {
                int64_t n = cur[w * W + d];
                counts[(static_cast<size_t>(p) * W + w) * W + d] = n;
                sum += n;
            }
            int64_t expect = 1;
            // binom(C,w) * 3^w
            int64_t b = 1;
            for (int i = 1; i <= w; ++i) b = b * (C - w + i) / i;
            expect = b * pow3_table[w];
            if (sum != expect) throw std::runtime_error("suffix orbit count mismatch");
        }
    }
    return counts;
}

inline int binom_cached(int n, int k, std::vector<int>& memo, int Nmax) {
    if (k < 0 || k > n) return 0;
    return memo[n * (Nmax + 1) + k];
}

inline KernelTable build_kernel(int N, int r, int m) {
    init_pows();
    KernelTable T;
    T.shape = make_shape(N, r, m);
    const int L = T.shape.L;
    const int C = T.shape.C;
    const int M = T.shape.m;
    if (L > 12) throw std::runtime_error("L too large for this builder");
    if (C > 20) throw std::runtime_error("C too large for this builder");
    T.n_pref = ipow_int(3, L);
    T.n_cent = T.n_pref * (C + 1);
    std::vector<int> pref(L);
    T.dmin.assign(T.n_pref, 0);
    std::vector<std::vector<int>> libs(r, std::vector<int>(L));
    for (int j = 0; j < r; ++j)
        for (int i = 0; i < L; ++i) libs[j][i] = (i < j ? 0 : 1);
    for (int id = 0; id < T.n_pref; ++id) {
        id_to_prefix(id, L, pref.data());
        int best = L + 1;
        for (int j = 0; j < r; ++j) {
            int d = 0;
            for (int i = 0; i < L; ++i) d += pref[i] != libs[j][i];
            if (d < best) best = d;
        }
        T.dmin[id] = static_cast<uint8_t>(best);
    }
    for (int id = 0; id < T.n_pref; ++id) {
        for (int w = 0; w <= C; ++w) {
            if (T.dmin[id] + w <= M) {
                T.point_pref.push_back(id);
                T.point_w.push_back(w);
            }
        }
    }
    T.n_pts = static_cast<int>(T.point_pref.size());
    const auto counts = suffix_counts(C);
    const int W = C + 1;
    std::vector<int> bin((C + 1) * (C + 1), 0);
    for (int n = 0; n <= C; ++n) {
        bin[n * W + 0] = 1;
        for (int k = 1; k <= n; ++k) {
            bin[n * W + k] = bin[(n - 1) * W + (k - 1)] + bin[(n - 1) * W + k];
        }
    }
    auto Cnk = [&](int n, int k) -> int64_t {
        if (k < 0 || k > n) return 0;
        return bin[n * W + k];
    };
    // prefix signature table: for every pair of prefix ids, forced and both
    std::vector<uint8_t> forced(static_cast<size_t>(T.n_pref) * T.n_pref);
    std::vector<uint8_t> both(static_cast<size_t>(T.n_pref) * T.n_pref);
    std::vector<int> A(L), B(L);
    for (int i = 0; i < T.n_pref; ++i) {
        id_to_prefix(i, L, A.data());
        for (int j = 0; j < T.n_pref; ++j) {
            id_to_prefix(j, L, B.data());
            int f = 0, bth = 0;
            for (int t = 0; t < L; ++t) {
                if (A[t] != B[t]) ++f;
                if (A[t] == 2 && B[t] == 2) ++bth;
            }
            forced[static_cast<size_t>(i) * T.n_pref + j] = static_cast<uint8_t>(f);
            both[static_cast<size_t>(i) * T.n_pref + j] = static_cast<uint8_t>(bth);
        }
    }
    // cumulative suffix counts
    std::vector<int64_t> cum(static_cast<size_t>(W) * W * (C + 2), 0);
    for (int p = 0; p <= C; ++p) {
        for (int w = 0; w <= C; ++w) {
            int64_t run = 0;
            cum[(static_cast<size_t>(p) * W + w) * (C + 2) + 0] = 0;
            for (int d = 0; d <= C; ++d) {
                run += counts[(static_cast<size_t>(p) * W + w) * W + d];
                cum[(static_cast<size_t>(p) * W + w) * (C + 2) + (d + 1)] = run;
            }
        }
    }
    T.rhs.assign(T.n_pts, 0);
    for (int pi = 0; pi < T.n_pts; ++pi) {
        int w = T.point_w[pi];
        // rhs = 2^L * binom(C,w) * 3^w
        __int128 v = static_cast<__int128>(pow2_table[L]) * Cnk(C, w) * pow3_table[w];
        if (v > static_cast<__int128>(INT64_MAX)) throw std::runtime_error("rhs overflow");
        T.rhs[pi] = static_cast<int64_t>(v);
    }
    T.scaled.assign(static_cast<size_t>(T.n_cent) * T.n_pts, 0);
    for (int pref = 0; pref < T.n_pref; ++pref) {
        for (int p = 0; p <= C; ++p) {
            int centre = pref * (C + 1) + p;
            for (int pi = 0; pi < T.n_pts; ++pi) {
                int qpref = T.point_pref[pi];
                int w = T.point_w[pi];
                int f = forced[static_cast<size_t>(pref) * T.n_pref + qpref];
                int bth = both[static_cast<size_t>(pref) * T.n_pref + qpref];
                int budget = M - f;
                if (budget < 0) continue;
                int64_t num = 0;
                for (int mis = 0; mis <= bth; ++mis) {
                    int end = budget - mis;
                    if (end < 0) continue;
                    if (end > C) end = C;
                    int64_t ways = Cnk(bth, mis);
                    int64_t got = cum[(static_cast<size_t>(p) * W + w) * (C + 2) + (end + 1)];
                    num += ways * got;
                }
                __int128 scaled = static_cast<__int128>(num) * pow2_table[L - bth];
                if (scaled > static_cast<__int128>(INT64_MAX)) throw std::runtime_error("scaled overflow");
                T.scaled[static_cast<size_t>(centre) * T.n_pts + pi] = static_cast<int64_t>(scaled);
            }
        }
    }
    return T;
}

inline std::string prefix_string(int id, int L) {
    int a[16];
    id_to_prefix(id, L, a);
    std::string s;
    for (int i = 0; i < L; ++i) s.push_back(a[i] == 2 ? 'X' : static_cast<char>('0' + a[i]));
    return s;
}
