///
///	\file   prime_moduli.cpp
///
///	\brief  Script file to generate and print the prime moduli
///			needed for the Ring LWE constructions. It takes in
///			the required number of bits k for the modulus and the
///			dimension n of the ring elements.
///
///			It outputs the largest (probable) prime number smaller
///			than 2^k which admits a negacyclic convolution for the
///			multiplication of elements of the ring Z_q / (x^n + 1).
///
///

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

constexpr int ITS = 1000000;

typedef __uint128_t uint128_t;

static inline uint128_t mulmod_u128(uint128_t a, uint128_t b, uint128_t m) {
  // Slow 128 bit multiplication: double-and-add without overflow
  // If need arises, use 256 bit package
  uint128_t res = 0;
  a %= m;
  while (b) {
    if (b & 1)
      res = (res >= m - a) ? res + a - m : res + a;
    b >>= 1;
    if (b)
      a = (a >= m - a) ? a + a - m : a + a;
  }
  return res;
}

static inline uint128_t powmod_u128(uint128_t a, uint128_t e, uint128_t m) {
  uint128_t r = 1;
  while (e) {
    if (e & 1)
      r = mulmod_u128(r, a, m);
    e >>= 1;
    if (e)
      a = mulmod_u128(a, a, m);
  }
  return r;
}

static bool miller_rabin_witness128(uint128_t a, uint128_t n) {
  if (a % n == 0)
    return true; // a ≡ 0 mod n --> “passes” this base
  // write n-1 = d * 2^s with d odd
  uint128_t d = n - 1, s = 0;
  while ((d & 1) == 0) {
    d >>= 1;
    ++s;
  }

  uint128_t x = powmod_u128(a, d, n);
  if (x == 1 || x == n - 1)
    return true;
  for (uint128_t r = 1; r < s; ++r) {
    x = mulmod_u128(x, x, n);
    if (x == n - 1)
      return true;
  }
  return false; // definitely composite for this base
}

bool is_prime_u128(uint128_t n) {
  if (n < 2)
    return false;
  // small primes quick path
  for (uint128_t p :
       {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull}) {
    if (n == p)
      return true;
    if (n % p == 0)
      return n == p;
  }
  // probabilistic MR for n < 2^128
  static constexpr std::array<uint64_t, 24> bases128 = {
      2ULL,  3ULL,  5ULL,  7ULL,  11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL,
      41ULL, 43ULL, 47ULL, 53ULL, 59ULL, 61ULL, 67ULL, 71ULL, 73ULL, 79ULL, 83ULL, 89ULL};
  for (uint128_t a : bases128) {
    if (!miller_rabin_witness128(a, n))
      return false;
  }
  return true;
}

// distinct prime factors of m
static inline std::vector<uint128_t> distinct_primes128(uint128_t m) {
  std::vector<uint128_t> p;
  if (m % 2 == 0) {
    p.push_back(2);
    while (m % 2 == 0)
      m /= 2;
  }
  for (uint128_t f = 3; f * f <= m; f += 2) {
    if (m % f == 0) {
      p.push_back(f);
      while (m % f == 0)
        m /= f;
    }
  }
  if (m > 1)
    p.push_back(m);
  return p;
}

// returns a primitive 2n-th root of unity mod prime q; assumes 2n | (q-1)
uint128_t primitive_root_2n128(uint128_t q, uint128_t n) {
  const uint128_t two_n = 2 * n;
  if ((q - 1) % two_n != 0)
    throw std::runtime_error("2n does not divide q-1");
  const uint128_t M = (q - 1) / two_n;
  const auto primes = distinct_primes128(two_n);

  std::mt19937_64 rng(std::random_device{});
  std::uniform_int_distribution<uint128_t> dist(2, q - 2);

  for (;;) {
    uint128_t a = dist(rng);
    uint128_t r = powmod_u128(a, M, q); // r in μ_{2n}
    if (r == 1)
      continue; // order 1 --> reject quickly
    bool ok = true;
    for (uint128_t p : primes) {
      if (powmod_u128(r, two_n / p, q) == 1) {
        ok = false;
        break;
      }
    }
    if (ok)
      return r; // r has exact order 2n
  }
}

//
// Printing functions
//

// Operator to print 128-bit numbers
std::ostream &operator<<(std::ostream &os, uint128_t value) {
  char buffer[129]; // It fits the entire value + the null terminated char
  char *d = std::end(buffer);

  // Put the null terminator
  --d;
  *d = '\0';

  do {
    --d;
    *d = "0123456789"[value % 10];
    value /= 10;
  } while (value != 0);

  std::string number(d);
  os << number;

  return os;
}

// Function to print 128-bit numbers in binary
std::string binStr(uint128_t input) {
  std::bitset<128> hi{static_cast<unsigned long long>(input >> 64)},
      lo{static_cast<unsigned long long>(input)}, bits{(hi << 64) | lo};

  return bits.to_string;
}

// Function to print a ruler above the bitset
std::string binRuler(string linePrefix) {
  std::ostringstream oss;

  // First line
  oss << linePrefix;
  int ctr = 127;
  while (ctr >= 0) {
    if (ctr % 10 == 9) {
      oss << "|";
    } else {
      oss << " ";
    };
    ctr--;
  }
  oss << endl;

  // Second line
  oss << linePrefix;
  ctr = 127;
  while (ctr >= 0) {
    oss << (ctr % 10);
    ctr--;
  }

  return oss.str;
}

// Function to print two 64 bit numbers in hex
std::string hexStr(uint128_t value) {
  std::ostringstream oss;

  char buffer[130]; // It fits the entire value + the extra chars
  char *d = std::end(buffer);

  // Put the null terminator
  --d;
  *d = '\0';

  int ctr = 64 / 4; // Split at 64 bits (each hex is 4 bits)
  do {
    --d;
    *d = "0123456789ABCDEF"[value % 16];
    if (--ctr == 0) {
      --d;
      *d = ' ';
    }
    value /= 16;
  } while (value != 0);

  std::string number(d);
  oss << number;

  return oss.str;
}

void PrintMaximumGoodPrime {
  uint64_t k;
  uint64_t n;
  cout << "=========================================" << endl;
  cout << "Give me the number of bits of the modulus k: ";
  cin >> k;
  cout << "Give me the dimension n: ";
  cin >> n;

  // Make sure than n is a power of 2
  if ((n & (n - 1)) != 0) {
    cout << "n is not a power of 2" << endl;
    return;
  }

  uint128_t qmax = 1;
  qmax <<= (uint128_t)k;
  qmax -= 1;

  cout << binRuler("                ") << endl;
  cout << "q can be up to " << "(" << binStr(qmax) << ")_2 = " << endl << qmax << endl;

  // q - 1 has to be divided by 2 * n
  uint128_t dn = 2 * n;
  uint128_t q = (qmax / dn) * dn + 1ull; // This is the maximum that is allowed

  for (int i = 0; i < ITS; i++) {
    if (is_prime_u128(q)) {
      cout << endl << "First probable prime: " << q << endl;
      cout << binRuler("") << endl;
      cout << binStr(q) << endl;
      break;
    } else {
      q -= dn;
    }
  }

  uint128_t psi = primitive_root_2n128(q, n);

  // Test
  uint128_t res = powmod_u128(psi, 2 * n, q);
  if (res == 1UL) {
    cout << endl << "Found root psi = " << psi << " res = " << res << endl;
  } else {
    cout << "Something went wrong!" << psi << " res = " << res << endl;
  }

  cout << "In Hex: q = " << hexStr(q) << endl << "      psi = " << hexStr(psi) << endl;
}

int main {

  PrintMaximumGoodPrime;

  return 0;
}
