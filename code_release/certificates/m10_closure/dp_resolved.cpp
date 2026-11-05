// Exact radius-m load of a resolved near-cut shell measure.
// Centres are ternary. Letter 2 stands for either non-binary letter.
// A shell is (anchor, guard symbols, resolved symbols, far edit counts).
// The reported maximum is the maximum of sum coeff * (shell words in the ball).

#include <algorithm>
#include <array>
#include <map>
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

struct Shell {
    int j;
    std::vector<int> gsym, resL, resR;
    int t, u, fl, fr;
    I coeff;
};

struct Layout {
    int g0, g1, bL, bR, start, end, farL, farR, wlen;
};

static Layout make_layout(int j, int N, int depth) {
    Layout L;
    L.g0 = std::max(0, j - 1);
    L.g1 = std::min(N, j + 1);
    int left_len = L.g0;
    int right_len = N - L.g1;
    L.bL = std::min(depth, left_len);
    L.bR = std::min(depth, right_len);
    L.start = L.g0 - L.bL;
    L.end = L.g1 + L.bR;
    L.farL = L.start;
    L.farR = N - L.end;
    L.wlen = L.end - L.start;
    return L;
}

static long long shell_size_of(const Shell& sh, int N, int m, int depth) {
    Layout L = make_layout(sh.j, N, depth);
    if ((int)sh.gsym.size() != L.g1 - L.g0) return 0;
    if ((int)sh.resL.size() != L.bL || (int)sh.resR.size() != L.bR) return 0;
    if (sh.t > L.farL || sh.u > L.farR || sh.fl > sh.t || sh.fr > sh.u) return 0;
    int edits = sh.t + sh.u;
    for (int s : sh.gsym) edits += s != 0;
    for (int s : sh.resL) edits += s != 0;
    for (int s : sh.resR) edits += s != 0;
    if (edits != m) return 0;
    long long ways = 1;
    for (int s : sh.gsym) if (s == 2) ways *= 2;
    for (int s : sh.resL) if (s == 2) ways *= 2;
    for (int s : sh.resR) if (s == 2) ways *= 2;
    auto choose = [](int n, int k) {
        if (k < 0 || k > n) return 0LL;
        long long r = 1;
        for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1);
        return r;
    };
    ways *= choose(L.farL, sh.t);
    ways *= choose(sh.t, sh.fl);
    ways *= choose(L.farR, sh.u);
    ways *= choose(sh.u, sh.fr);
    ways *= 1LL << (sh.t - sh.fl + sh.u - sh.fr);
    return ways;
}

int main(int argc, char** argv) {
    try {
        if (argc < 3) return 2;
        std::string mode = argv[1];
        std::ifstream in(argv[2]);
        int N, m, depth, flip_cap, free_guard, nshell;
        in >> N >> m >> depth >> flip_cap >> free_guard >> nshell;
        if (!in || N < 1 || N > 26 || m < 0 || m > 12 || depth < 0 || depth > 4 ||
            flip_cap < 0 || flip_cap > m || nshell < 0 || nshell > 200000)
            return 3;
        std::vector<Shell> shells;
        shells.reserve(nshell);
        for (int i = 0; i < nshell; i++) {
            Shell sh;
            int nG, nL, nR;
            in >> sh.j >> nG >> nL >> nR;
            sh.gsym.resize(nG);
            sh.resL.resize(nL);
            sh.resR.resize(nR);
            for (int k = 0; k < nG; k++) in >> sh.gsym[k];
            for (int k = 0; k < nL; k++) in >> sh.resL[k];
            for (int k = 0; k < nR; k++) in >> sh.resR[k];
            std::string cs;
            in >> sh.t >> sh.u >> sh.fl >> sh.fr >> cs;
            if (!in) return 4;
            sh.coeff = parse_int(cs);
            Layout L = make_layout(sh.j, N, depth);
            if (sh.j < 0 || sh.j > N || nG != L.g1 - L.g0 || nL != L.bL || nR != L.bR) return 4;
            if (!free_guard) {
                for (int s : sh.gsym)
                    if (s != 0) return 4;
            }
            shells.push_back(std::move(sh));
        }

        const int W = 2 * depth + 2;
        const int Q = m + 1;
        const int FMAX = flip_cap;
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

        if (mode == "count") {
            int ncent;
            in >> ncent;
            std::vector<std::vector<int>> centres(ncent, std::vector<int>(N));
            for (int i = 0; i < ncent; i++) {
                std::string s;
                in >> s;
                if ((int)s.size() != N) return 5;
                for (int k = 0; k < N; k++) centres[i][k] = s[k] - '0';
            }
            for (int si = 0; si < nshell; si++) {
                const Shell& sh = shells[si];
                Layout L = make_layout(sh.j, N, depth);
                long long sz = shell_size_of(sh, N, m, depth);
                for (int ci = 0; ci < ncent; ci++) {
                    const auto& c = centres[ci];
                    int cur[13] = {0};
                    cur[0] = 1;
                    auto convolve = [&](auto site_cost) {
                        int nxt[13] = {0};
                        for (int d = 0; d <= m; d++) if (cur[d]) {
                            for (int e = 0; e <= m - d; e++) if (site_cost[e])
                                nxt[d + e] += cur[d] * site_cost[e];
                        }
                        for (int d = 0; d <= m; d++) cur[d] = nxt[d];
                    };
                    auto cost_of = [&](int symbol, int letter, int pat) {
                        int w[13] = {0};
                        if (symbol == 0) w[letter == pat ? 0 : 1] = 1;
                        else if (symbol == 1) w[letter == (1 - pat) ? 0 : 1] = 1;
                        else if (letter >= 2) { w[0] = 1; w[1] = 1; }
                        else w[1] = 2;
                        convolve(w);
                    };
                    int gi = 0;
                    for (int pos = L.g0; pos < L.g1; pos++, gi++)
                        cost_of(sh.gsym[gi], c[pos], pos < sh.j ? 1 : 0);
                    for (int k = 0; k < L.bL; k++) {
                        int pos = L.g0 - 1 - k;
                        cost_of(sh.resL[k], c[pos], 1);
                    }
                    for (int k = 0; k < L.bR; k++) {
                        int pos = L.g1 + k;
                        cost_of(sh.resR[k], c[pos], 0);
                    }
                    auto arm = [&](int from, int to, int agree, int edits, int flips) -> std::vector<long long> {
                        int ag = 0, dg = 0, cg = 0;
                        for (int p = from; p < to; p++) {
                            if (c[p] == agree) ag++;
                            else if (c[p] < 2) dg++;
                            else cg++;
                        }
                        const auto& row = poly[pid(ag, dg, cg)];
                        std::vector<long long> pref(Q + 1, 0);
                        if (row.empty()) return pref;
                        for (int d = 0; d <= m; d++) pref[d + 1] = pref[d] + row[pidx(edits, flips, d)];
                        return pref;
                    };
                    auto left = arm(0, L.start, 1, sh.t, sh.fl);
                    auto right = arm(L.end, N, 0, sh.u, sh.fr);
                    long long total = 0;
                    for (int d = 0; d <= m; d++) if (cur[d]) {
                        int limit = m - d;
                        for (int dl = 0; dl <= limit; dl++) {
                            long long lv = left[dl + 1] - left[dl];
                            if (lv) total += cur[d] * lv * right[limit - dl + 1];
                        }
                    }
                    std::cout << total << " " << sz << "\n";
                }
            }
            return 0;
        }

        if (mode == "matrix") {
            if (argc != 4) return 2;
            int ncent;
            in >> ncent;
            if (!in || ncent < 1 || ncent > 20000) return 5;
            std::vector<std::vector<int>> centres(ncent, std::vector<int>(N));
            for (int i = 0; i < ncent; i++) {
                std::string s;
                in >> s;
                if ((int)s.size() != N) return 5;
                for (int k = 0; k < N; k++) centres[i][k] = s[k] - '0';
            }
            std::vector<double> mat((size_t)ncent * nshell, 0.0);
            for (int si = 0; si < nshell; si++) {
                const Shell& sh = shells[si];
                Layout L = make_layout(sh.j, N, depth);
                long long sz = shell_size_of(sh, N, m, depth);
                if (sz <= 0) return 4;
                for (int ci = 0; ci < ncent; ci++) {
                    const auto& c = centres[ci];
                    int cur[13] = {0};
                    cur[0] = 1;
                    auto convolve = [&](int symbol, int letter, int pat) {
                        int wgt[13] = {0};
                        if (symbol == 0) wgt[letter == pat ? 0 : 1] = 1;
                        else if (symbol == 1) wgt[letter == (1 - pat) ? 0 : 1] = 1;
                        else if (letter >= 2) { wgt[0] = 1; wgt[1] = 1; }
                        else wgt[1] = 2;
                        int nxt[13] = {0};
                        for (int d = 0; d <= m; d++) if (cur[d])
                            for (int e = 0; e <= m - d; e++) if (wgt[e])
                                nxt[d + e] += cur[d] * wgt[e];
                        for (int d = 0; d <= m; d++) cur[d] = nxt[d];
                    };
                    int gi = 0;
                    for (int pos = L.g0; pos < L.g1; pos++, gi++)
                        convolve(sh.gsym[gi], c[pos], pos < sh.j ? 1 : 0);
                    for (int k = 0; k < L.bL; k++)
                        convolve(sh.resL[k], c[L.g0 - 1 - k], 1);
                    for (int k = 0; k < L.bR; k++)
                        convolve(sh.resR[k], c[L.g1 + k], 0);
                    auto arm = [&](int from, int to, int agree, int edits, int flips, long long* pref) {
                        int ag = 0, dg = 0, cg = 0;
                        for (int p = from; p < to; p++) {
                            if (c[p] == agree) ag++;
                            else if (c[p] < 2) dg++;
                            else cg++;
                        }
                        const auto& row = poly[pid(ag, dg, cg)];
                        for (int d = 0; d <= m; d++) pref[d + 1] = pref[d];
                        if (row.empty()) return;
                        for (int d = 0; d <= m; d++) pref[d + 1] = pref[d] + row[pidx(edits, flips, d)];
                    };
                    long long left[14] = {0}, right[14] = {0};
                    arm(0, L.start, 1, sh.t, sh.fl, left);
                    arm(L.end, N, 0, sh.u, sh.fr, right);
                    long long total = 0;
                    for (int d = 0; d <= m; d++) if (cur[d]) {
                        int limit = m - d;
                        for (int dl = 0; dl <= limit; dl++) {
                            long long lv = left[dl + 1] - left[dl];
                            if (lv) total += cur[d] * lv * right[limit - dl + 1];
                        }
                    }
                    mat[(size_t)ci * nshell + si] = (double)total / (double)sz;
                }
            }
            std::ofstream out(argv[3], std::ios::binary);
            out.write(reinterpret_cast<const char*>(mat.data()), (std::streamsize)(mat.size() * sizeof(double)));
            if (!out) return 6;
            std::cout << "matrix " << ncent << " " << nshell << "\n";
            return 0;
        }

        if (mode != "max") return 2;

        const int HIST = 1;
        std::vector<int> pow3(W + 1, 1);
        for (int k = 1; k <= W; k++) pow3[k] = pow3[k - 1] * 3;
        const int HCAP = pow3[W];

        struct AnchorPack {
            Layout lay;
            int nL, nR;
            std::vector<I> cost;
        };
        std::vector<AnchorPack> packs(N + 1);
        for (int j = 0; j <= N; j++) {
            Layout L = make_layout(j, N, depth);
            AnchorPack ap;
            ap.lay = L;
            ap.nL = L.farL + 1;
            ap.nR = L.farR + 1;
            ap.cost.assign((size_t)pow3[L.wlen] * ap.nL * ap.nL * ap.nR * ap.nR, 0);
            packs[j] = std::move(ap);
        }

        auto left_id = [](int ag, int dg, int nL) { return ag * nL + dg; };

        if (std::getenv("DP_REFERENCE_PRECOMPUTE")) {
        for (const Shell& sh : shells) {
            if (sh.coeff == 0) continue;
            Layout L = make_layout(sh.j, N, depth);
            long long sz = shell_size_of(sh, N, m, depth);
            if (sz <= 0) throw std::runtime_error("empty shell");
            AnchorPack& ap = packs[sh.j];
            int nL = ap.nL, nR = ap.nR;
            int wlen = L.wlen;
            int nW = pow3[wlen];
            // Aggregate this shell into the anchor table.
            // Window id enumerates centre letters on [start, end).
            for (int wid = 0; wid < nW; wid++) {
                int letters[16];
                int tmp = wid;
                for (int k = 0; k < wlen; k++) {
                    letters[k] = tmp % 3;
                    tmp /= 3;
                }
                // digits were stored least-significant first, so letters[0] is the
                // leftmost site only if we encoded it that way. Encode wid with
                // the leftmost site as the least significant digit: letters[0] is start.
                int cur[13] = {0};
                cur[0] = 1;
                auto convolve = [&](int symbol, int letter, int pat) {
                    int wgt[13] = {0};
                    if (symbol == 0) wgt[letter == pat ? 0 : 1] = 1;
                    else if (symbol == 1) wgt[letter == (1 - pat) ? 0 : 1] = 1;
                    else if (letter >= 2) { wgt[0] = 1; wgt[1] = 1; }
                    else wgt[1] = 2;
                    int nxt[13] = {0};
                    for (int d = 0; d <= m; d++) if (cur[d])
                        for (int e = 0; e <= m - d; e++) if (wgt[e])
                            nxt[d + e] += cur[d] * wgt[e];
                    for (int d = 0; d <= m; d++) cur[d] = nxt[d];
                };
                for (int k = 0; k < L.g1 - L.g0; k++) {
                    int pos = L.g0 + k;
                    int letter = letters[pos - L.start];
                    convolve(sh.gsym[k], letter, pos < sh.j ? 1 : 0);
                }
                for (int k = 0; k < L.bL; k++) {
                    int pos = L.g0 - 1 - k;
                    convolve(sh.resL[k], letters[pos - L.start], 1);
                }
                for (int k = 0; k < L.bR; k++) {
                    int pos = L.g1 + k;
                    convolve(sh.resR[k], letters[pos - L.start], 0);
                }
                for (int ag = 0; ag <= L.farL; ag++) {
                    for (int dg = 0; dg <= L.farL - ag; dg++) {
                        int cg = L.farL - ag - dg;
                        const auto& Lp = poly[pid(ag, dg, cg)];
                        if (Lp.empty()) continue;
                        for (int ar = 0; ar <= L.farR; ar++) {
                            for (int dr = 0; dr <= L.farR - ar; dr++) {
                                int cr = L.farR - ar - dr;
                                const auto& Rp = poly[pid(ar, dr, cr)];
                                if (Rp.empty()) continue;
                                long long ways = 0;
                                for (int d0 = 0; d0 <= m; d0++) if (cur[d0]) {
                                    int limit = m - d0;
                                    for (int dl = 0; dl <= limit; dl++) {
                                        long long lv = Lp[pidx(sh.t, sh.fl, dl)];
                                        if (!lv) continue;
                                        long long rv = 0;
                                        // sum of right poly over drd <= limit-dl
                                        for (int drd = 0; drd <= limit - dl; drd++)
                                            rv += Rp[pidx(sh.u, sh.fr, drd)];
                                        ways += cur[d0] * lv * rv;
                                    }
                                }
                                if (!ways) continue;
                                size_t idx = ((size_t)wid * nL * nL + (size_t)left_id(ag, dg, nL)) * nR * nR
                                    + (size_t)left_id(ar, dr, nR);
                                ap.cost[idx] += sh.coeff * (I)ways;
                            }
                        }
                    }
                }
            }
        }


        } else {

        // At a fixed anchor, far-arm counts depend only on (t,u,fl,fr).
        // Sum coefficient-weighted window polynomials before their convolution
        // with those far-arm counts. All summands and operations are integers.
        using Key = std::array<int,5>;
        std::map<Key, std::vector<const Shell*>> groups;
        for (const Shell& sh : shells) {
            if (!sh.coeff) continue;
            if (shell_size_of(sh, N, m, depth) <= 0)
                throw std::runtime_error("empty shell");
            groups[{sh.j,sh.t,sh.u,sh.fl,sh.fr}].push_back(&sh);
        }
        for (const auto& entry : groups) {
            const Key& key = entry.first;
            int j=key[0], t=key[1], u=key[2], fl=key[3], fr=key[4];
            Layout L=make_layout(j,N,depth);
            AnchorPack& ap=packs[j];
            int nL=ap.nL,nR=ap.nR,nW=pow3[L.wlen];
            std::vector<I> local((size_t)nW*Q,0);
            for (const Shell* sp : entry.second) {
                const Shell& sh=*sp;
                for (int wid=0;wid<nW;wid++) {
                    int letters[16],tmp=wid;
                    for (int k=0;k<L.wlen;k++) {letters[k]=tmp%3;tmp/=3;}
                    int cur[13]={0};cur[0]=1;
                    auto convolve=[&](int symbol,int letter,int pat) {
                        int weights[13]={0},next[13]={0};
                        if (symbol==0) weights[letter==pat?0:1]=1;
                        else if (symbol==1) weights[letter==(1-pat)?0:1]=1;
                        else if (letter>=2) {weights[0]=1;weights[1]=1;}
                        else weights[1]=2;
                        for (int d=0;d<=m;d++) if(cur[d])
                            for(int e=0;e<=m-d;e++) if(weights[e])
                                next[d+e]+=cur[d]*weights[e];
                        std::copy(next,next+13,cur);
                    };
                    for(int k=0;k<L.g1-L.g0;k++) {
                        int pos=L.g0+k;
                        convolve(sh.gsym[k],letters[pos-L.start],pos<j?1:0);
                    }
                    for(int k=0;k<L.bL;k++) {
                        int pos=L.g0-1-k;
                        convolve(sh.resL[k],letters[pos-L.start],1);
                    }
                    for(int k=0;k<L.bR;k++) {
                        int pos=L.g1+k;
                        convolve(sh.resR[k],letters[pos-L.start],0);
                    }
                    for(int d=0;d<=m;d++)
                        local[(size_t)wid*Q+d]+=sh.coeff*(I)cur[d];
                }
            }
            for(int ag=0;ag<=L.farL;ag++) for(int dg=0;dg<=L.farL-ag;dg++) {
                const auto& lp=poly[pid(ag,dg,L.farL-ag-dg)];
                if(lp.empty()) continue;
                for(int ar=0;ar<=L.farR;ar++) for(int dr=0;dr<=L.farR-ar;dr++) {
                    const auto& rp=poly[pid(ar,dr,L.farR-ar-dr)];
                    if(rp.empty()) continue;
                    I prefix[13]={0};
                    for(int dl=0;dl<=m;dl++) {
                        long long lv=lp[pidx(t,fl,dl)];
                        if(!lv) continue;
                        for(int rd=0;rd<=m-dl;rd++)
                            prefix[dl+rd]+=(I)lv*(I)rp[pidx(u,fr,rd)];
                    }
                    for(int d=1;d<=m;d++) prefix[d]+=prefix[d-1];
                    for(int wid=0;wid<nW;wid++) {
                        I ways=0;
                        for(int d=0;d<=m;d++)
                            ways+=local[(size_t)wid*Q+d]*prefix[m-d];
                        size_t index=((size_t)wid*nL*nL+(size_t)left_id(ag,dg,nL))*nR*nR
                                      +(size_t)left_id(ar,dr,nR);
                        ap.cost[index]+=ways;
                    }
                }
            }
        }

        }

        // Sweep. Reuse one anchor's table.
        std::vector<U> p3(N + 1, 1);
        for (int k = 1; k <= N; k++) p3[k] = p3[k - 1] * 3ULL;

        struct State { I value; U word; bool ok; };
        const int nZO = (N + 1) * (N + 2) / 2;
        std::vector<I> best_num(nZO, -1);
        std::vector<U> best_word(nZO, 0);
        std::vector<std::pair<int, int>> tasks;
        for (int Z = 0; Z <= N; Z++)
            for (int O = 0; O <= N - Z; O++) tasks.push_back({Z, O});

        unsigned nthreads = 6;
        if (const char* env = std::getenv("DP_THREADS")) {
            int t = std::atoi(env);
            if (t >= 1 && t <= 6) nthreads = (unsigned)t;
        }
        if (nthreads > tasks.size()) nthreads = (unsigned)tasks.size();
        std::atomic<int> cursor{0};
        auto worker = [&]() {
            std::vector<State> prev((size_t)SIDE * SIDE * HCAP), next(prev.size());
            auto sid = [&](int z, int o, int h) { return (z * SIDE + o) * HCAP + h; };
            while (true) {
                int job = cursor.fetch_add(1);
                if (job >= (int)tasks.size()) break;
                int Z = tasks[job].first;
                int O = tasks[job].second;
                int X = N - Z - O;
                for (auto& st : prev) st.ok = false;
                prev[sid(0, 0, 0)] = State{0, 0, true};
                for (int p = 0; p < N; p++) {
                    for (auto& st : next) st.ok = false;
                    for (int z = 0; z <= Z && z <= p; z++) {
                        for (int o = 0; o <= O && z + o <= p; o++) {
                            int x = p - z - o;
                            if (x > X || x < 0) continue;
                            for (int h = 0; h < HCAP; h++) {
                                State st = prev[sid(z, o, h)];
                                if (!st.ok) continue;
                                for (int cn = 0; cn < 3; cn++) {
                                    int zz = z + (cn == 0);
                                    int oo = o + (cn == 1);
                                    int xx = x + (cn == 2);
                                    if (zz > Z || oo > O || xx > X) continue;
                                    U word = st.word * 3ULL + (U)cn;
                                    I val = st.value;
                                    int np = p + 1;
                                    // Score every anchor whose window ends at np.
                                    for (int j = 0; j <= N; j++) {
                                        const Layout& L = packs[j].lay;
                                        if (L.end != np) continue;
                                        int wlen = L.wlen;
                                        int wid = 0, mul = 1;
                                        int wz = 0, wo = 0, wx = 0;
                                        for (int k = 0; k < wlen; k++) {
                                            int letter = (int)((word / p3[np - 1 - (L.start + k)]) % 3ULL);
                                            wid += letter * mul;
                                            mul *= 3;
                                            wz += letter == 0;
                                            wo += letter == 1;
                                            wx += letter == 2;
                                        }
                                        int ag = oo - wo;
                                        int dg = zz - wz;
                                        int cg = (np - zz - oo) - wx;
                                        int ar = Z - zz;
                                        int dr = O - oo;
                                        int cr = X - (np - zz - oo);
                                        if (ag < 0 || dg < 0 || cg < 0 || ar < 0 || dr < 0 || cr < 0)
                                            throw std::runtime_error("negative arm count");
                                        if (ag + dg + cg != L.farL || ar + dr + cr != L.farR)
                                            throw std::runtime_error("arm length");
                                        int nL = packs[j].nL, nR = packs[j].nR;
                                        size_t idx = ((size_t)wid * nL * nL + (size_t)left_id(ag, dg, nL)) * nR * nR
                                            + (size_t)left_id(ar, dr, nR);
                                        val += packs[j].cost[idx];
                                    }
                                    int nh = (int)(word % (U)HCAP);
                                    // For p+1 < W the modulus still equals the word, since word < 3^{p+1} <= 3^W.
                                    State& dst = next[sid(zz, oo, nh)];
                                    if (!dst.ok || val > dst.value || (val == dst.value && word < dst.word))
                                        dst = State{val, word, true};
                                }
                            }
                        }
                    }
                    prev.swap(next);
                }
                I best = -1;
                U bw = 0;
                bool any = false;
                for (int h = 0; h < HCAP; h++) {
                    State st = prev[sid(Z, O, h)];
                    if (!st.ok) continue;
                    any = true;
                    if (st.value > best || (st.value == best && st.word < bw)) {
                        best = st.value;
                        bw = st.word;
                    }
                }
                if (!any) throw std::runtime_error("missing class");
                int slot = 0;
                for (int z = 0; z < Z; z++) slot += N - z + 1;
                slot += O;
                best_num[slot] = best;
                best_word[slot] = bw;
            }
        };
        auto t0 = std::chrono::steady_clock::now();
        std::vector<std::thread> pool;
        for (unsigned t = 0; t < nthreads; t++) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
        I global = -1;
        U gw = 0;
        int slot = 0;
        for (int Z = 0; Z <= N; Z++) {
            for (int O = 0; O <= N - Z; O++, slot++) {
                I best = best_num[slot];
                U bw = best_word[slot];
                if (std::getenv("DP_PRINT_TOTALS"))
                    std::cout << "total " << Z << " " << O << " " << dec(best) << " " << bw << "\n";
                if (best > global || (best == global && bw < gw)) {
                    global = best;
                    gw = bw;
                }
            }
        }
        std::string letters(N, '0');
        U tmp = gw;
        for (int k = N - 1; k >= 0; k--) {
            letters[k] = char('0' + int(tmp % 3ULL));
            tmp /= 3ULL;
        }
        std::cout << "maximum " << dec(global) << "\n";
        std::cout << "witness " << letters << "\n";
        std::cout << "classes " << dec((I)p3[N]) << "\n";
        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "seconds " << sec << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 6;
    }
}
