// Exact maximum radius-m load of a tagged guarded-shell measure.
// Centres are ternary words over {0,1,2}. Letter 2 stands for either
// non-binary letter. Each shell coefficient is an integer. The reported
// maximum is the maximum of sum coeff * (shell words inside the ball).
//
// A shell is (g, j, t, u, fl, fr, tagL, tagR). The guard
// [j-g, j+g) matches the pattern. tagL and tagR are the cells immediately
// outside the guard, or -1 if that cell does not exist. Tag 0 is fixed to
// the pattern, tag 1 is the opposite binary letter, tag 2 is a non-binary
// letter. The far arms carry t and u edits, of which fl and fr are flips.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using U = uint64_t;

struct U256 {
    U w[4]{};

    bool is_zero() const { return (w[0] | w[1] | w[2] | w[3]) == 0; }

    static U256 parse(const std::string& s) {
        if (s.empty()) throw std::runtime_error("empty integer");
        U256 a;
        for (char c : s) {
            if (c < '0' || c > '9') throw std::runtime_error("nonnegative integer required");
            unsigned __int128 carry = static_cast<unsigned>(c - '0');
            for (int i = 0; i < 4; i++) {
                unsigned __int128 cur = static_cast<unsigned __int128>(a.w[i]) * 10 + carry;
                a.w[i] = static_cast<U>(cur);
                carry = cur >> 64;
            }
            if (carry) throw std::runtime_error("coefficient exceeds 256 bits");
        }
        return a;
    }

    void add_mul(const U256& coeff, U bucket) {
        if (bucket == 0 || coeff.is_zero()) return;
        unsigned __int128 carry = 0;
        for (int i = 0; i < 4; i++) {
            unsigned __int128 cur = static_cast<unsigned __int128>(w[i])
                + static_cast<unsigned __int128>(coeff.w[i]) * bucket
                + carry;
            w[i] = static_cast<U>(cur);
            carry = cur >> 64;
        }
        if (carry) throw std::runtime_error("load exceeds 256 bits");
    }

    void add(const U256& other) {
        unsigned __int128 carry = 0;
        for (int i = 0; i < 4; i++) {
            unsigned __int128 cur = static_cast<unsigned __int128>(w[i])
                + other.w[i] + carry;
            w[i] = static_cast<U>(cur);
            carry = cur >> 64;
        }
        if (carry) throw std::runtime_error("load exceeds 256 bits");
    }

    int cmp(const U256& o) const {
        for (int i = 3; i >= 0; i--) {
            if (w[i] < o.w[i]) return -1;
            if (w[i] > o.w[i]) return 1;
        }
        return 0;
    }
};

static std::string dec(U256 a) {
    if (a.is_zero()) return "0";
    std::string s;
    while (!a.is_zero()) {
        unsigned __int128 rem = 0;
        for (int i = 3; i >= 0; i--) {
            unsigned __int128 cur = (rem << 64) | a.w[i];
            a.w[i] = static_cast<U>(cur / 10);
            rem = cur % 10;
        }
        s.push_back(static_cast<char>('0' + static_cast<int>(rem)));
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct State {
    U256 value;
    U word = 0;
    bool has = false;
};

struct Shell {
    int g, j, t, u, fl, fr, tagL, tagR;
    U256 coeff;
};

struct Family {
    int g = 0;
    bool hasL = false;
    bool hasR = false;
    int nTL = 1;
    int nTR = 1;
    std::vector<int> left_len, right_len, nE;
    std::vector<std::vector<U256>> cost;
};

struct Opt {
    int dist;
    int ways;
};

static int tag_options(int tag, int centre, int pattern_letter, Opt* out) {
    if (tag < 0) {
        out[0] = {0, 1};
        return 1;
    }
    if (tag == 0) {
        out[0] = {centre == pattern_letter ? 0 : 1, 1};
        return 1;
    }
    if (tag == 1) {
        int flip = 1 - pattern_letter;
        out[0] = {centre == flip ? 0 : 1, 1};
        return 1;
    }
    if (centre >= 2) {
        out[0] = {0, 1};
        out[1] = {1, 1};
        return 2;
    }
    out[0] = {1, 2};
    return 1;
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) return 2;
        std::ifstream in(argv[1]);
        int N, m, nshell;
        in >> N >> m >> nshell;
        if (!in || N < 1 || N > 24 || m < 0 || m > 12 || nshell < 1 || nshell > 100000) return 3;
        std::vector<Shell> shells;
        shells.reserve(nshell);
        int G = 1;
        int FMAX = 0;
        for (int i = 0; i < nshell; i++) {
            int g, j, t, u, fl, fr, tagL, tagR;
            std::string cs;
            in >> g >> j >> t >> u >> fl >> fr >> tagL >> tagR >> cs;
            if (!in || g < 1 || g > 4 || j < 0 || j > N || t < 0 || u < 0 || fl < 0 || fr < 0)
                return 4;
            if (fl > t || fr > u || tagL < -1 || tagL > 2 || tagR < -1 || tagR > 2) return 4;
            int extra = (tagL == 1 || tagL == 2) + (tagR == 1 || tagR == 2);
            if (t + u + extra != m) return 4;
            bool hasL = j - g - 1 >= 0;
            bool hasR = j + g < N;
            if ((tagL >= 0) != hasL || (tagR >= 0) != hasR) return 4;
            G = std::max(G, g);
            FMAX = std::max(FMAX, std::max(fl, fr));
            shells.push_back({g, j, t, u, fl, fr, tagL, tagR, U256::parse(cs)});
        }
        const int Q = m + 1;
        const int SIDE = N + 1;
        auto pid = [&](int a, int b, int c) { return (a * SIDE + b) * SIDE + c; };
        auto pidx = [&](int s, int f, int d) { return (s * (FMAX + 1) + f) * Q + d; };
        const int PLEN = Q * (FMAX + 1) * Q;
        std::vector<std::vector<long long>> poly(static_cast<size_t>(SIDE) * SIDE * SIDE);
        poly[pid(0, 0, 0)].assign(PLEN, 0);
        poly[pid(0, 0, 0)][pidx(0, 0, 0)] = 1;
        auto add_to = [&](std::vector<long long>& dst, int s, int f, int d, long long v) {
            if (s < 0 || f < 0 || d < 0 || s > m || f > FMAX || d > m || v == 0) return;
            dst[pidx(s, f, d)] += v;
        };
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

        std::vector<Family> fams;
        // One family object per guard width. Cost tables are per boundary.
        {
            std::vector<int> gs;
            for (auto& sh : shells) {
                if (std::find(gs.begin(), gs.end(), sh.g) == gs.end()) gs.push_back(sh.g);
            }
            std::sort(gs.begin(), gs.end());
            for (int g : gs) {
                Family fam;
                fam.g = g;
                fam.left_len.assign(N + 1, 0);
                fam.right_len.assign(N + 1, 0);
                fam.nE.assign(N + 1, 0);
                fam.cost.assign(N + 1, {});
                for (int j = 0; j <= N; j++) {
                    bool hasL = j - g - 1 >= 0;
                    bool hasR = j + g < N;
                    int L = hasL ? (j - g - 1) : std::max(0, j - g);
                    int R = hasR ? (N - (j + g + 1)) : std::max(0, N - (j + g));
                    fam.left_len[j] = L;
                    fam.right_len[j] = R;
                    int nE = 2 * g + 1;
                    fam.nE[j] = nE;
                    int nTL = hasL ? 3 : 1;
                    int nTR = hasR ? 3 : 1;
                    size_t cells = static_cast<size_t>(L + 1) * (L + 1) * (R + 1) * (R + 1) * nE * nTL * nTR;
                    fam.cost[j].assign(cells, U256{});
                }
                // Geometry is the same for every boundary of this guard, but
                // hasL/hasR depend on j. The lookup below recomputes them.
                fams.push_back(std::move(fam));
            }
        }
        auto family = [&](int g) -> Family& {
            for (auto& f : fams) if (f.g == g) return f;
            throw std::runtime_error("missing family");
        };

        for (auto& sh : shells) {
            if (sh.coeff.is_zero()) continue;
            Family& fam = family(sh.g);
            int j = sh.j;
            int L = fam.left_len[j];
            int R = fam.right_len[j];
            int nE = fam.nE[j];
            bool hasL = j - sh.g - 1 >= 0;
            bool hasR = j + sh.g < N;
            int nTL = hasL ? 3 : 1;
            int nTR = hasR ? 3 : 1;
            int nRtypes = (R + 1) * (R + 1);
            int nLspan = L + 1;
            int nRspan = R + 1;
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
                            unsigned __int128 conv_wide[16];
                            std::memset(conv_wide, 0, sizeof(conv_wide));
                            bool any = false;
                            for (int dl = 0; dl <= m; dl++) {
                                long long lv = Lp[pidx(sh.t, sh.fl, dl)];
                                if (!lv) continue;
                                for (int drd = 0; drd + dl <= m; drd++) {
                                    long long rv = Rp[pidx(sh.u, sh.fr, drd)];
                                    if (!rv) continue;
                                    conv_wide[dl + drd] += static_cast<unsigned __int128>(lv) * static_cast<unsigned __int128>(rv);
                                    any = true;
                                }
                            }
                            if (!any) continue;
                            U conv[16];
                            for (int s = 0; s <= m; s++) {
                                if (conv_wide[s] > static_cast<unsigned __int128>(~U(0)))
                                    throw std::runtime_error("arm product exceeds 64 bits");
                                conv[s] = static_cast<U>(conv_wide[s]);
                            }
                            int cLmax = hasL ? 3 : 1;
                            int cRmax = hasR ? 3 : 1;
                            for (int cL = 0; cL < cLmax; cL++) {
                                for (int cR = 0; cR < cRmax; cR++) {
                                    Opt lo[2], ro[2];
                                    int nl = tag_options(sh.tagL, hasL ? cL : 0, 1, lo);
                                    int nr = tag_options(sh.tagR, hasR ? cR : 0, 0, ro);
                                    U temp[16];
                                    std::memset(temp, 0, sizeof(temp));
                                    for (int il = 0; il < nl; il++) {
                                        for (int ir = 0; ir < nr; ir++) {
                                            int extra = lo[il].dist + ro[ir].dist;
                                            unsigned __int128 ways = static_cast<unsigned __int128>(lo[il].ways) * ro[ir].ways;
                                            for (int s = 0; s + extra <= m; s++) {
                                                if (!conv[s]) continue;
                                                unsigned __int128 piece = conv[s] * ways;
                                                for (int e = 0; e < nE && s + extra + e <= m; e++) {
                                                    unsigned __int128 sum = temp[e] + piece;
                                                    if (sum > static_cast<unsigned __int128>(~U(0)))
                                                        throw std::runtime_error("bucket exceeds 64 bits");
                                                    temp[e] = static_cast<U>(sum);
                                                }
                                            }
                                        }
                                    }
                                    size_t type_base = (static_cast<size_t>(ag * nLspan + dg) * nRtypes
                                        + static_cast<size_t>(ar * nRspan + dr)) * nE;
                                    for (int e = 0; e < nE; e++) {
                                        if (!temp[e]) continue;
                                        size_t at = ((type_base + e) * nTL + cL) * nTR + cR;
                                        fam.cost[j][at].add_mul(sh.coeff, temp[e]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        std::vector<U> pow3(N + 1, 1);
        for (int k = 1; k <= N; k++) pow3[k] = pow3[k - 1] * 3;
        U hist_base = 1;
        for (int k = 0; k < 2 * G + 1; k++) hist_base *= 3;
        if (hist_base > 200000) throw std::runtime_error("history too large");
        const int HIST = static_cast<int>(hist_base);
        auto letter_at = [&](U word, int length, int k) -> int {
            if (k < 0 || k >= length) return -1;
            return static_cast<int>((word / pow3[length - 1 - k]) % 3);
        };
        auto term = [&](int j, int length, int z, int o, U word, int Z, int O) -> U256 {
            if (j < 0 || j > N || j > length) throw std::runtime_error("bad boundary");
            int pref[3] = {z, o, length - z - o};
            int seen[3] = {0, 0, 0};
            for (int k = j; k < length; k++) {
                int a = letter_at(word, length, k);
                if (a < 0) throw std::runtime_error("missing letter");
                seen[a]++;
            }
            int head[3] = {pref[0] - seen[0], pref[1] - seen[1], pref[2] - seen[2]};
            int X = N - Z - O;
            U256 value;
            for (auto& fam : fams) {
                int g = fam.g;
                int left[3] = {head[0], head[1], head[2]};
                int right[3] = {Z - head[0], O - head[1], X - head[2]};
                int e = 0;
                int L0 = std::max(0, j - g);
                for (int k = L0; k < j; k++) {
                    int a = letter_at(word, length, k);
                    if (a < 0) throw std::runtime_error("guard letter");
                    left[a]--;
                    e += (a != 1);
                }
                int R1 = std::min(N, j + g);
                for (int k = j; k < R1; k++) {
                    int a = letter_at(word, length, k);
                    if (a < 0) throw std::runtime_error("guard letter");
                    right[a]--;
                    e += (a != 0);
                }
                bool hasL = j - g - 1 >= 0;
                bool hasR = j + g < N;
                int cL = 0, cR = 0;
                if (hasL) {
                    int k = j - g - 1;
                    int a = letter_at(word, length, k);
                    if (a < 0) throw std::runtime_error("left tag");
                    left[a]--;
                    cL = a;
                }
                if (hasR) {
                    int k = j + g;
                    int a = letter_at(word, length, k);
                    if (a < 0) throw std::runtime_error("right tag");
                    right[a]--;
                    cR = a;
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
                int nTL = hasL ? 3 : 1;
                int nTR = hasR ? 3 : 1;
                int nE = fam.nE[j];
                int nRtypes = (R + 1) * (R + 1);
                size_t at = ((static_cast<size_t>(ag * (L + 1) + dg) * nRtypes
                    + static_cast<size_t>(ar * (R + 1) + dr)) * nE + e) * nTL + cL;
                at = at * nTR + cR;
                value.add(fam.cost[j][at]);
            }
            return value;
        };

        const int nZO = (N + 1) * (N + 2) / 2;
        std::vector<U256> best_num(nZO);
        std::vector<U> best_word(nZO, 0);
        std::vector<char> best_has(nZO, 0);
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
        std::atomic<int> cursor{0};
        std::atomic<bool> bad{false};
        std::string error;
        std::mutex error_mu;
        auto fail = [&](const std::string& msg) {
            bool expected = false;
            if (bad.compare_exchange_strong(expected, true)) {
                std::lock_guard<std::mutex> lock(error_mu);
                error = msg;
            }
        };
        std::vector<std::thread> pool;
        auto worker = [&]() {
            try {
                std::vector<State> prev(static_cast<size_t>(SIDE) * SIDE * HIST), next(prev.size());
                auto sid = [&](int z, int o, int h) { return (z * SIDE + o) * HIST + h; };
                while (!bad.load()) {
                    int job = cursor.fetch_add(1);
                    if (job >= static_cast<int>(tasks.size())) break;
                    int Z = tasks[job].first;
                    int O = tasks[job].second;
                    int X = N - Z - O;
                    std::fill(prev.begin(), prev.end(), State{});
                    prev[sid(0, 0, 0)] = State{U256{}, 0, true};
                    for (int p = 0; p < N; p++) {
                        std::fill(next.begin(), next.end(), State{});
                        for (int z = 0; z <= Z && z <= p; z++) {
                            for (int o = 0; o <= O && z + o <= p; o++) {
                                int x = p - z - o;
                                if (x > X) continue;
                                for (int h = 0; h < HIST; h++) {
                                    State st = prev[sid(z, o, h)];
                                    if (!st.has) continue;
                                    for (int cn = 0; cn < 3; cn++) {
                                        int zz = z + (cn == 0);
                                        int oo = o + (cn == 1);
                                        int xx = x + (cn == 2);
                                        if (zz > Z || oo > O || xx > X) continue;
                                        U word = st.word * 3 + static_cast<U>(cn);
                                        U256 val = st.value;
                                        int j = p - G;
                                        if (j >= 0) val.add(term(j, p + 1, zz, oo, word, Z, O));
                                        int nh = static_cast<int>(word % static_cast<U>(HIST));
                                        State& dst = next[sid(zz, oo, nh)];
                                        if (!dst.has || val.cmp(dst.value) > 0 || (val.cmp(dst.value) == 0 && word < dst.word))
                                            dst = State{val, word, true};
                                    }
                                }
                            }
                        }
                        prev.swap(next);
                    }
                    U256 best;
                    U bw = 0;
                    bool has = false;
                    for (int h = 0; h < HIST; h++) {
                        State st = prev[sid(Z, O, h)];
                        if (!st.has) continue;
                        U256 val = st.value;
                        for (int j = N - G; j <= N; j++) val.add(term(j, N, Z, O, st.word, Z, O));
                        if (!has || val.cmp(best) > 0 || (val.cmp(best) == 0 && st.word < bw)) {
                            best = val;
                            bw = st.word;
                            has = true;
                        }
                    }
                    if (!has) throw std::runtime_error("missing class");
                    int slot = 0;
                    for (int z = 0; z < Z; z++) slot += N - z + 1;
                    slot += O;
                    best_num[slot] = best;
                    best_word[slot] = bw;
                    best_has[slot] = 1;
                }
            } catch (const std::exception& ex) {
                fail(ex.what());
            }
        };
        auto t0 = std::chrono::steady_clock::now();
        for (unsigned t = 0; t < nthreads; t++) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
        if (bad.load()) {
            std::cerr << error << "\n";
            return 6;
        }
        U256 global;
        U gw = 0;
        bool ghas = false;
        int slot = 0;
        for (int Z = 0; Z <= N; Z++) {
            for (int O = 0; O <= N - Z; O++, slot++) {
                if (!best_has[slot]) throw std::runtime_error("missing class output");
                U256 best = best_num[slot];
                U bw = best_word[slot];
                std::cout << "total " << Z << " " << O << " " << dec(best) << " " << bw << "\n";
                if (!ghas || best.cmp(global) > 0 || (best.cmp(global) == 0 && bw < gw)) {
                    global = best;
                    gw = bw;
                    ghas = true;
                }
            }
        }
        std::string letters(N, '0');
        U tmp = gw;
        for (int k = N - 1; k >= 0; k--) {
            letters[k] = static_cast<char>('0' + int(tmp % 3));
            tmp /= 3;
        }
        std::cout << "maximum " << dec(global) << " " << gw << "\n";
        std::cout << "witness " << letters << "\n";
        std::cout << "classes " << pow3[N] << "\n";
        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "seconds " << sec << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 6;
    }
}
