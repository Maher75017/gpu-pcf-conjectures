// cf.cu — recherche massive de fractions continues polynomiales pour des constantes difficiles.
// V = a0 + b1/(a1 + b2/(a2 + ...)),  a_n = A3 n^3 + A2 n^2 + A1 n + A0,  b_n = c * forme(n).
// Identification : V = (p + q C)/(r + s C), |p,q,r,s| <= 12.
// Contrôles positifs dans l'espace : Apéry 6/zeta(3) (34,51,27,5 ; -n^6) et 5/zeta(2) (0,11,11,3 ; +n^4).
// Contrôle négatif : 5 constantes aléatoires dans la même table.
// GPU : nvcc -O3 -arch=native -o cf cf.cu && ./cf        CPU (test) : g++ -O2 -x c++ -o cf_cpu cf.cu && ./cf_cpu
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif
HD inline bool finite_(double x) { return x == x && fabs(x) <= 1.7976931348623157e308; }

#define NTERMS 200
#ifndef A3N
#define A3N 81   // A3 in 0..80
#endif
#ifndef A2N
#define A2N 129  // A2 in 0..128
#endif
#ifndef A1N
#define A1N 81   // A1 in 0..80
#endif
#ifndef A0N
#define A0N 32   // A0 in 1..32
#endif
#define NFORM 8
#define NCS 32   // c in -16..-1, 1..16
#define K 12
#define TOL 3e-12
#define MAXH 2000000

struct Entry { double v; int32_t id; };
struct Hit { uint64_t idx; int32_t id; double V; };

#define NCONST 25
static const char *CNAME[NCONST] = {"zeta3", "zeta5", "zeta7", "catalan", "gamma", "zeta2", "pi", "e", "ln2", "ln3", "ln5",
  "pi3", "pi_ln2", "ln2_2", "pi_sqrt3", "pi_sqrt2", "ln_1sqrt2", "lemniscate", "gauss", "pi_e",
  "alea0", "alea1", "alea2", "alea3", "alea4"};
static double CVAL[NCONST];
static void init_consts() {
  double pi = 3.14159265358979323846, l2 = log(2.0), a = 1, b = sqrt(2.0);
  for (int i = 0; i < 40; i++) { double t = (a + b) / 2; b = sqrt(a * b); a = t; }
  double v[NCONST] = {1.2020569031595942854, 1.0369277551433699263, 1.0083492773819228268, 0.91596559417721901505,
    0.57721566490153286061, pi * pi / 6, pi, exp(1.0), l2, log(3.0), log(5.0), pi * pi * pi, pi * l2, l2 * l2,
    pi / sqrt(3.0), pi / sqrt(2.0), log(1 + sqrt(2.0)), tgamma(0.25) * tgamma(0.25) / (2 * sqrt(2 * pi)), 1 / a, pi * exp(1.0),
    1.3846920711942651, 0.7263159348230112, 2.1730581239740027, 1.0571937621503329, 0.6189041732014821};
  for (int i = 0; i < NCONST; i++) CVAL[i] = v[i];
}

HD inline void decode(uint64_t i, int &A3, int &A2, int &A1, int &A0, int &f, int &c) {
  int ci = i % NCS; i /= NCS; c = ci < 16 ? ci - 16 : ci - 15;
  f = i % NFORM; i /= NFORM;
  A0 = 1 + (int)(i % A0N); i /= A0N;
  A1 = i % A1N; i /= A1N;
  A2 = i % A2N; i /= A2N;
  A3 = (int)i;
}
HD inline double bval(int f, int c, double n) {
  double t = 2 * n - 1, n2 = n * n;
  switch (f) {
    case 0: return c * n2 * n2 * n2;          // c n^6
    case 1: return c * n2 * n2;               // c n^4
    case 2: return c * n2 * t * t;            // c n^2 (2n-1)^2
    case 3: return c * n2 * n2 * t * t;       // c n^4 (2n-1)^2
    case 4: return c * t * t * t * t;         // c (2n-1)^4
    case 5: return c * n2 * n * t * t * t;    // c n^3 (2n-1)^3
    case 6: return c * t * t * t * t * t * t; // c (2n-1)^6
    default: return c * n2 * n2 * (n + 1) * (n + 1); // c n^4 (n+1)^2
  }
}
// renvoie true si convergée ; V en sortie
HD inline bool evalcf(int A3, int A2, int A1, int A0, int f, int c, double &V) {
  double p0 = 1, p1 = A0, q0 = 0, q1 = 1, vm = 0;
  for (int k = 1; k <= NTERMS; k++) {
    double n = k, a = ((A3 * n + A2) * n + A1) * n + A0, b = bval(f, c, n);
    double p2 = a * p1 + b * p0, q2 = a * q1 + b * q0;
    p0 = p1; p1 = p2; q0 = q1; q1 = q2;
    double m = fmax(fabs(p1), fabs(q1));
    if (!(m > 0) || !finite_(m)) return false;
    if (m > 1e100) { double s = 1.0 / m; p0 *= s; p1 *= s; q0 *= s; q1 *= s; }
    if (k == NTERMS / 2) { if (q1 == 0) return false; vm = p1 / q1; }
  }
  if (q1 == 0) return false;
  V = p1 / q1;
  if (!finite_(V) || fabs(V) < 1e-6 || fabs(V) > 1e6) return false;
  return fabs(V - vm) <= 1e-13 * fabs(V);
}
HD inline void match(const Entry *T, int nT, uint64_t idx, double V, Hit *H, unsigned long long *nh) {
  double lo = V - TOL * fabs(V), hi = V + TOL * fabs(V);
  int a = 0, b = nT;
  while (a < b) { int m = (a + b) / 2; if (T[m].v < lo) a = m + 1; else b = m; }
  for (int j = a; j < nT && T[j].v <= hi; j++) {
#ifdef __CUDA_ARCH__
    unsigned long long w = atomicAdd(nh, 1ULL);
#else
    unsigned long long w = (*nh)++;
#endif
    if (w < MAXH) { H[w].idx = idx; H[w].id = T[j].id; H[w].V = V; }
  }
}

static int gcd4(int a, int b, int c, int d) { return std::gcd(std::gcd(abs(a), abs(b)), std::gcd(abs(c), abs(d))); }
static std::vector<Entry> build_table() {
  std::vector<Entry> T;
  init_consts();
  for (int ci = 0; ci < NCONST; ci++)
    for (int p = -K; p <= K; p++) for (int q = -K; q <= K; q++)
      for (int r = -K; r <= K; r++) for (int s = -K; s <= K; s++) {
        if (q * r - p * s == 0 || (q == 0 && s == 0) || gcd4(p, q, r, s) != 1) continue;
        if (r < 0 || (r == 0 && s < 0)) continue;           // signe normalisé
        double C = CVAL[ci], den = r + s * C;
        if (fabs(den) < 1e-9) continue;
        double v = (p + q * C) / den;
        if (fabs(v) < 1e-6 || fabs(v) > 1e6) continue;
        int id = ci * 390625 + (((p + K) * 25 + (q + K)) * 25 + (r + K)) * 25 + (s + K);
        T.push_back({v, id});
      }
  std::sort(T.begin(), T.end(), [](const Entry &x, const Entry &y) { return x.v < y.v; });
  return T;
}
static void write_hits(const Hit *H, unsigned long long nh, const char *fn) {
  FILE *f = fopen(fn, "w");
  unsigned long long n = nh < MAXH ? nh : MAXH;
  for (unsigned long long i = 0; i < n; i++) {
    int A3, A2, A1, A0, fo, c; decode(H[i].idx, A3, A2, A1, A0, fo, c);
    int id = H[i].id, ci = id / 390625, rest = id % 390625;
    int s = rest % 25 - K; rest /= 25; int r = rest % 25 - K; rest /= 25; int q = rest % 25 - K; int p = rest / 25 - K;
    fprintf(f, "%d %d %d %d %d %d %s %d %d %d %d %.17g\n", A3, A2, A1, A0, fo, c, CNAME[ci], p, q, r, s, H[i].V);
  }
  fclose(f);
  printf("correspondances : %llu (ecrites dans %s)\n", nh, fn);
}
#define TOTAL ((uint64_t)A3N * A2N * A1N * A0N * NFORM * NCS)

#ifdef __CUDACC__
__global__ void kern(const Entry *T, int nT, Hit *H, unsigned long long *nh, unsigned long long *nconv, uint64_t i0, uint64_t i1) {
  unsigned long long conv = 0;
  for (uint64_t i = i0 + blockIdx.x * (uint64_t)blockDim.x + threadIdx.x; i < i1; i += (uint64_t)gridDim.x * blockDim.x) {
    int A3, A2, A1, A0, f, c; decode(i, A3, A2, A1, A0, f, c);
    double V;
    if (evalcf(A3, A2, A1, A0, f, c, V)) { conv++; match(T, nT, i, V, H, nh); }
  }
  atomicAdd(nconv, conv);
}
int main() {
  auto T = build_table();
  printf("espace : %llu fractions continues ; table : %zu valeurs\n", (unsigned long long)TOTAL, T.size());
  Entry *dT; Hit *dH; unsigned long long *dn, *dc, zero = 0;
  cudaMalloc(&dT, T.size() * sizeof(Entry)); cudaMalloc(&dH, MAXH * sizeof(Hit));
  cudaMalloc(&dn, 8); cudaMalloc(&dc, 8);
  cudaMemcpy(dT, T.data(), T.size() * sizeof(Entry), cudaMemcpyHostToDevice);
  cudaMemcpy(dn, &zero, 8, cudaMemcpyHostToDevice); cudaMemcpy(dc, &zero, 8, cudaMemcpyHostToDevice);
  unsigned long long nh, nc;
  cudaEvent_t e0, e1; cudaEventCreate(&e0); cudaEventCreate(&e1); cudaEventRecord(e0);
  const uint64_t CH = 250000000ULL;
  for (uint64_t i0 = 0; i0 < TOTAL; i0 += CH) {
    uint64_t i1 = i0 + CH < TOTAL ? i0 + CH : TOTAL;
    kern<<<170 * 16, 256>>>(dT, (int)T.size(), dH, dn, dc, i0, i1);
    cudaDeviceSynchronize(); cudaMemcpy(&nh, dn, 8, cudaMemcpyDeviceToHost);
    printf("  %5.1f %%  correspondances : %llu\n", 100.0 * i1 / TOTAL, nh); fflush(stdout);
  }
  cudaEventRecord(e1); cudaEventSynchronize(e1); float ms; cudaEventElapsedTime(&ms, e0, e1);
  cudaError_t err = cudaGetLastError(); if (err != cudaSuccess) { printf("erreur CUDA : %s\n", cudaGetErrorString(err)); return 1; }
  cudaMemcpy(&nh, dn, 8, cudaMemcpyDeviceToHost); cudaMemcpy(&nc, dc, 8, cudaMemcpyDeviceToHost);
  std::vector<Hit> H(MAXH); cudaMemcpy(H.data(), dH, MAXH * sizeof(Hit), cudaMemcpyDeviceToHost);
  printf("convergees : %llu ; temps GPU : %.1f s\n", nc, ms / 1000);
  write_hits(H.data(), nh, "hits.txt");
  return 0;
}
#else
int main(int argc, char **argv) {  // test CPU : contrôles + échantillon
  auto T = build_table();
  printf("table : %zu valeurs\n", T.size());
  std::vector<Hit> H(MAXH); unsigned long long nh = 0, nc = 0;
  auto idx_of = [](int A3, int A2, int A1, int A0, int f, int c) {
    uint64_t i = A3; i = i * A2N + A2; i = i * A1N + A1; i = i * A0N + (A0 - 1); i = i * NFORM + f;
    return i * NCS + (c < 0 ? c + 16 : c + 15); };
  std::vector<uint64_t> L = {idx_of(34, 51, 27, 5, 0, -1), idx_of(0, 11, 11, 3, 1, 1)};
  uint64_t step = TOTAL / (argc > 1 ? atoll(argv[1]) : 200000);
  for (uint64_t i = 7; i < TOTAL; i += step) L.push_back(i);
  for (uint64_t i : L) { int A3, A2, A1, A0, f, c; decode(i, A3, A2, A1, A0, f, c); double V;
    if (evalcf(A3, A2, A1, A0, f, c, V)) { nc++; match(T.data(), (int)T.size(), i, V, H.data(), &nh); } }
  printf("testees : %zu ; convergees : %llu\n", L.size(), nc);
  write_hits(H.data(), nh, "hits_cpu.txt");
  return 0;
}
#endif
