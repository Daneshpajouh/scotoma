// Exact maximum radius-m load of a guarded-shell measure.
// Centres are ternary words over {0,1,2}. Letter 2 stands for either
// non-binary letter. The caller passes integer shell coefficients
// coeff = v * (H / S). The reported maximum is the maximum of
// sum coeff * (number of shell words inside the ball).

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using I = __int128_t;
using U = uint64_t;

static std::string dec(I x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x) {
        s.push_back(char('0' + int(x % 10)));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

static I parse_int(const std::string& s) {
    if (s.empty()) throw std::runtime_error("empty integer");
    I v = 0;
    for (char c : s) {
        if (c < '0' || c > '9') throw std::runtime_error("nonnegative integer required");
        v = v * 10 + (c - '0');
    }
    return v;
}

struct State {
    I value = -1;
    U word = 0;
};

struct Shell {
    int g, j, t, fl, fr;
    I coeff;
};

struct Family {
    int g;
    int flips;
    std::vector<int> left_len, right_len;
    std::vector<std::vector<I>> cost;
    std::vector<int> stride_r, stride_e, nL, nR;
};

int main(int argc, char** argv) {
    try {
        if (argc != 2) return 2;
        std::ifstream in(argv[1]);
        int N, m, nshell;
        in >> N >> m >> nshell;
        if (!in || N < 1 || N > 24 || m < 0 || m > 12 || nshell < 0 || nshell > 20000) return 3;
        std::vector<Shell> shells;
        shells.reserve(nshell);
        int G = 1;
        for (int i = 0; i < nshell; i++) {
            int g, j, t, fl, fr;
            std::string cs;
            in >> g >> j >> t >> fl >> fr >> cs;
            if (!in || g < 1 || g > 6 || j < 0 || j > N || t < 0 || t > m || fl < 0 || fr < 0)
                return 4;
            G = std::max(G, g);
            shells.push_back({g, j, t, fl, fr, parse_int(cs)});
        }
        const int FMAX = G;  // flip index 0..G-1 is enough; shells use <= g-1
        const int Q = m + 1;
        const int SIDE = N + 1;
        auto pid = [&](int a, int b, int c) { return (a * SIDE + b) * SIDE + c; };
        auto pidx = [&](int s, int f, int d) { return (s * (FMAX + 1) + f) * Q + d; };
        const int PLEN = Q * (FMAX + 1) * Q;
        std::vector<std::vector<long long>> poly((SIDE) * (SIDE) * (SIDE));
        poly[pid(0, 0, 0)].assign(PLEN, 0);
        poly[pid(0, 0, 0)][pidx(0, 0, 0)] = 1;
        auto add_to = [&](std::vector<long long>& dst, int s, int f, int d, long long v) {
            if (s < 0 || f < 0 || d < 0 || s > m || f > FMAX || d > m || v == 0) return;
            dst[pidx(s, f, d)] += v;
        };
        // One canonical parent per composition. Summing over every order of
        // position types would multiply the polynomial by a multinomial.
        for (int len = 1; len <= N; len++) {
            for (int a = 0; a <= len; a++) {
                for (int b = 0; b <= len - a; b++) {
                    int c = len - a - b;
                    int kind, pa, pb, pc;
                    if (c > 0) {
                        kind = 2;
                        pa = a;
                        pb = b;
                        pc = c - 1;
                    } else if (b > 0) {
                        kind = 1;
                        pa = a;
                        pb = b - 1;
                        pc = c;
                    } else {
                        kind = 0;
                        pa = a - 1;
                        pb = b;
                        pc = c;
                    }
                    const auto& src = poly[pid(pa, pb, pc)];
                    if (src.empty()) throw std::runtime_error("missing arm parent");
                    auto& dst = poly[pid(a, b, c)];
                    dst.assign(PLEN, 0);
                    for (int s = 0; s <= m; s++) {
                        for (int f = 0; f <= FMAX; f++) {
                            for (int d = 0; d <= m; d++) {
                                long long v = src[pidx(s, f, d)];
                                if (!v) continue;
                                if (kind == 0) {
                                    add_to(dst, s, f, d, v);
                                    add_to(dst, s + 1, f, d + 1, 2 * v);
                                    add_to(dst, s + 1, f + 1, d + 1, v);
                                } else if (kind == 1) {
                                    add_to(dst, s, f, d + 1, v);
                                    add_to(dst, s + 1, f, d + 1, 2 * v);
                                    add_to(dst, s + 1, f + 1, d, v);
                                } else {
                                    add_to(dst, s, f, d + 1, v);
                                    add_to(dst, s + 1, f, d, v);
                                    add_to(dst, s + 1, f, d + 1, v);
                                    add_to(dst, s + 1, f + 1, d + 1, v);
                                }
                            }
                        }
                    }
                }
            }
        }
        std::vector<int> gs;
        for (auto& sh : shells) {
            if (std::find(gs.begin(), gs.end(), sh.g) == gs.end()) gs.push_back(sh.g);
        }
        std::sort(gs.begin(), gs.end());
        std::vector<Family> fams;
        for (int g : gs) {
            Family fam;
            fam.g = g;
            fam.flips = g - 1;
            fam.left_len.assign(N + 1, 0);
            fam.right_len.assign(N + 1, 0);
            fam.cost.assign(N + 1, {});
            fam.nL.assign(N + 1, 0);
            fam.nR.assign(N + 1, 0);
            fam.stride_r.assign(N + 1, 0);
            fam.stride_e.assign(N + 1, 0);
            for (int j = 0; j <= N; j++) {
                int L = std::max(0, j - g);
                int R = std::max(0, N - j - g);
                fam.left_len[j] = L;
                fam.right_len[j] = R;
                int nL = (L + 1) * (L + 1);
                int nR = (R + 1) * (R + 1);
                int nE = 2 * g + 1;
                fam.nL[j] = nL;
                fam.nR[j] = nR;
                fam.stride_e[j] = 1;
                fam.stride_r[j] = nE;
                fam.cost[j].assign(static_cast<size_t>(nL) * nR * nE, 0);
            }
            fams.push_back(std::move(fam));
        }
        auto fam_of = [&](int g) -> Family& {
            for (auto& f : fams) if (f.g == g) return f;
            throw std::runtime_error("missing family");
        };
        auto ab_index = [](int a, int b, int len) {
            return a * (len + 1) + b;
        };
        for (auto& sh : shells) {
            if (sh.coeff == 0) continue;
            Family& fam = fam_of(sh.g);
            int j = sh.j;
            int L = fam.left_len[j];
            int R = fam.right_len[j];
            int nE = 2 * sh.g + 1;
            int nR = fam.nR[j];
            for (int ag = 0; ag <= L; ag++) {
                for (int dg = 0; dg <= L - ag; dg++) {
                    int cg = L - ag - dg;
                    const auto& Lp = poly[pid(ag, dg, cg)];
                    if (Lp.empty()) continue;
                    for (int ar = 0; ar <= R; ar++) {
                        for (int dr = 0; dr <= R - ar; dr++) {
                            int cr = R - ar - dr;
                            const auto& Rp = poly[pid(ar, dr, cr)];
                            if (Rp.empty()) continue;
                            long long bucket[16];
                            std::memset(bucket, 0, sizeof(bucket));
                            for (int dl = 0; dl <= m; dl++) {
                                long long lv = Lp[pidx(sh.t, sh.fl, dl)];
                                if (!lv) continue;
                                for (int drd = 0; drd + dl <= m; drd++) {
                                    long long rv = Rp[pidx(m - sh.t, sh.fr, drd)];
                                    if (!rv) continue;
                                    long long prod = lv * rv;
                                    for (int e = 0; e < nE && dl + drd + e <= m; e++) bucket[e] += prod;
                                }
                            }
                            size_t base = (static_cast<size_t>(ab_index(ag, dg, L)) * nR + ab_index(ar, dr, R)) * nE;
                            for (int e = 0; e < nE; e++) {
                                if (bucket[e]) fam.cost[j][base + e] += sh.coeff * static_cast<I>(bucket[e]);
                            }
                        }
                    }
                }
            }
        }

        std::vector<U> pow3(N + 1, 1);
        for (int k = 1; k <= N; k++) pow3[k] = pow3[k - 1] * 3;
        const int HIST = static_cast<int>(pow3[std::max(0, 2 * G - 1)]);
        auto letter_at = [&](U word, int p, int k) -> int {
            if (k < 0 || k >= p) return -1;
            return static_cast<int>((word / pow3[p - 1 - k]) % 3);
        };
        auto term = [&](int j, int p, int z, int o, U word, int Z, int O) -> I {
            if (j < 0 || j > N || j > p) throw std::runtime_error("bad boundary");
            int pref[3] = {z, o, p - z - o};
            int seen[3] = {0, 0, 0};
            for (int k = j; k < p; k++) {
                int a = letter_at(word, p, k);
                if (a < 0) throw std::runtime_error("missing letter");
                seen[a]++;
            }
            int head[3] = {pref[0] - seen[0], pref[1] - seen[1], pref[2] - seen[2]};
            int X = N - Z - O;
            I value = 0;
            for (auto& fam : fams) {
                int g = fam.g;
                int left[3] = {head[0], head[1], head[2]};
                int right[3] = {Z - head[0], O - head[1], X - head[2]};
                int e = 0;
                int L0 = std::max(0, j - g);
                for (int k = L0; k < j; k++) {
                    int a = letter_at(word, p, k);
                    if (a < 0) throw std::runtime_error("guard letter");
                    left[a]--;
                    e += (a != 1);
                }
                int R1 = std::min(N, j + g);
                for (int k = j; k < R1; k++) {
                    int a = letter_at(word, p, k);
                    if (a < 0) throw std::runtime_error("guard letter");
                    right[a]--;
                    e += (a != 0);
                }
                int L = fam.left_len[j];
                int R = fam.right_len[j];
                int ag = left[1], dg = left[0];
                int ar = right[0], dr = right[1];
                if (ag < 0 || dg < 0 || left[2] < 0 || ar < 0 || dr < 0 || right[2] < 0)
                    throw std::runtime_error("negative arm");
                if (ag + dg + left[2] != L || ar + dr + right[2] != R)
                    throw std::runtime_error("arm length");
                if (e < 0 || e > 2 * g) throw std::runtime_error("guard distance");
                int nE = 2 * g + 1;
                int nR = fam.nR[j];
                size_t base = (static_cast<size_t>(ag * (L + 1) + dg) * nR + (ar * (R + 1) + dr)) * nE;
                value += fam.cost[j][base + e];
            }
            return value;
        };

        const int nZO = (N + 1) * (N + 2) / 2;
        std::vector<I> best_num(nZO, -1);
        std::vector<U> best_word(nZO, 0);
        std::vector<std::pair<int, int>> tasks;
        tasks.reserve(nZO);
        for (int Z = 0; Z <= N; Z++)
            for (int O = 0; O <= N - Z; O++) tasks.push_back({Z, O});

        unsigned nthreads = std::thread::hardware_concurrency();
        if (nthreads == 0 || nthreads > 8) nthreads = 8;
        if (const char* env = std::getenv("DP_THREADS")) {
            int t = std::atoi(env);
            if (t >= 1 && t <= 8) nthreads = static_cast<unsigned>(t);
        }
        if (nthreads > tasks.size()) nthreads = static_cast<unsigned>(tasks.size());
        std::vector<std::thread> pool;
        std::atomic<int> cursor{0};
        auto worker = [&]() {
            std::vector<State> prev(static_cast<size_t>(SIDE) * SIDE * HIST), next(prev.size());
            auto sid = [&](int z, int o, int h) { return (z * SIDE + o) * HIST + h; };
            while (true) {
                int job = cursor.fetch_add(1);
                if (job >= static_cast<int>(tasks.size())) break;
                int Z = tasks[job].first;
                int O = tasks[job].second;
                int X = N - Z - O;
                std::fill(prev.begin(), prev.end(), State{});
                prev[sid(0, 0, 0)] = State{0, 0};
                for (int p = 0; p < N; p++) {
                    std::fill(next.begin(), next.end(), State{});
                    for (int z = 0; z <= Z && z <= p; z++) {
                        for (int o = 0; o <= O && z + o <= p; o++) {
                            int x = p - z - o;
                            if (x > X) continue;
                            for (int h = 0; h < HIST; h++) {
                                State st = prev[sid(z, o, h)];
                                if (st.value < 0) continue;
                                for (int cn = 0; cn < 3; cn++) {
                                    int zz = z + (cn == 0);
                                    int oo = o + (cn == 1);
                                    int xx = x + (cn == 2);
                                    if (zz > Z || oo > O || xx > X) continue;
                                    U word = st.word * 3 + static_cast<U>(cn);
                                    I val = st.value;
                                    int j = p + 1 - G;
                                    if (j >= 0) val += term(j, p + 1, zz, oo, word, Z, O);
                                    int nh = static_cast<int>(word % static_cast<U>(HIST));
                                    State& dst = next[sid(zz, oo, nh)];
                                    if (val > dst.value || (val == dst.value && word < dst.word))
                                        dst = State{val, word};
                                }
                            }
                        }
                    }
                    prev.swap(next);
                }
                I best = -1;
                U bw = 0;
                for (int h = 0; h < HIST; h++) {
                    State st = prev[sid(Z, O, h)];
                    if (st.value < 0) continue;
                    I val = st.value;
                    for (int j = N - G + 1; j <= N; j++) val += term(j, N, Z, O, st.word, Z, O);
                    if (val > best || (val == best && st.word < bw)) {
                        best = val;
                        bw = st.word;
                    }
                }
                if (best < 0) throw std::runtime_error("missing class");
                int slot = 0;
                for (int z = 0; z < Z; z++) slot += N - z + 1;
                slot += O;
                best_num[slot] = best;
                best_word[slot] = bw;
            }
        };
        auto t0 = std::chrono::steady_clock::now();
        for (unsigned t = 0; t < nthreads; t++) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
        I global = -1;
        U gw = 0;
        int slot = 0;
        for (int Z = 0; Z <= N; Z++) {
            for (int O = 0; O <= N - Z; O++, slot++) {
                I best = best_num[slot];
                U bw = best_word[slot];
                std::cout << "total " << Z << " " << O << " " << dec(best) << " " << bw << "\n";
                if (best > global || (best == global && bw < gw)) {
                    global = best;
                    gw = bw;
                }
            }
        }
        std::string letters;
        U tmp = gw;
        letters.assign(N, '0');
        for (int k = N - 1; k >= 0; k--) {
            letters[k] = char('0' + int(tmp % 3));
            tmp /= 3;
        }
        std::cout << "maximum " << dec(global) << " " << gw << "\n";
        std::cout << "witness " << letters << "\n";
        std::cout << "classes " << dec(static_cast<I>(pow3[N])) << "\n";
        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "seconds " << sec << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 6;
    }
}
