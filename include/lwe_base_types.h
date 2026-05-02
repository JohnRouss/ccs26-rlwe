///
///	\file   lwe_base_types.h
///
///	\brief  Header file to specify base types for Ring LWE
///         constructions. These constructions use rings
///         isomorphic to Z_q[x] / (x^n + 1) and fields
///         isomorphic to R[x] / (x^n + 1) with n a power
///         of 2.
///         We specify the base types for the ring and field
///         coefficients and the indexes they might admit.
///

// Comment below to switch to 128-bit Rcf's
#define __RCF64__

/// @brief Type of ring coefficients
#ifdef __RCF64__
typedef int64_t Rcf;
typedef uint64_t URcf;
#else
typedef __int128_t Rcf;
typedef __uint128_t URcf;
#endif

static constexpr URcf URCF_MAX = URcf(Rcf(-1L));
static constexpr Rcf RCF_MAX = URCF_MAX >> 1;

/// @brief Type of the dimension n for the rings, and the
///        logarithms of various quantities
typedef uint32_t Dim;

/// @brief Type of field coefficients and the complex ones
///        for their spectrums
typedef double Fcf;
typedef std::complex<Fcf> FScf;

// Functions to print and get the decimal length of an Rcf
#ifndef __RCF64__
std::ostream &operator<<(std::ostream &os, Rcf value);
#endif
size_t decimalLength(Rcf value);

// -------------------------------------------------------

// This function takes in any integer and returns the same
// modulo q in the range [-(q-1)/2, (q-1)/2].
static inline Rcf RRange(Rcf a, Rcf q) {
  a %= q;                    // a --> [-(q-1), (q-1)]
  a += q * (Rcf)(a < 0);     // In case a < 0 --> a \in [0,q-1]
  a -= q * (Rcf)(a > q / 2); // In case a > q/2 --> a \in [-(q-1)/2, (q-1)/2]
  return a;
}

// Modular operations over Z_q. They take in and return results in
// the range [-(q-1)/2, (q-1)/2] (q is an (positive) odd prime).
static inline Rcf RAdd(Rcf a, Rcf b, Rcf q) {
  // (a + b) does not overflow as each operand is in [-(q-1)/2, (q-1)/2].
  // Thus the sum is in the range [-(q-1), (q-1)].
  return RRange(a + b, q);
}

static inline Rcf RSub(Rcf a, Rcf b, Rcf q) {
  // (a - b) does not overflow as each operand is in [-(q-1)/2, (q-1)/2].
  // Thus the sum is in the range [-(q-1), (q-1)].
  return RRange(a - b, q);
}

static inline Rcf RMul(Rcf a, Rcf b, Rcf q) {
#ifdef __RCF64__
  __int128_t res = ((__int128_t)a * (__int128_t)b) % (__int128_t)q;
  return RRange((Rcf)res, q);
#else
  // Here a*b can overflow. In order to avoid overflows
  // we add to the result a, 2a, 4a,... (mod q) depending
  // on the bits of b.
  // For the algorithm to be constant time, we always
  // continue the process for all bits of q.

  // Get a positive b --> [q-(q-1)/2, q+(q-1)/2]
  b += q * (Rcf)(b < 0);
  URcf ub = (URcf)b;

  Rcf res = 0;
  Rcf tmp = 0;
  Rcf qCtr = q;

  // To achieve constant time (k steps where
  // k is the number of bits of q), we shift
  // q one bit to the right every iteration.
  while (qCtr > 0) {
    tmp = a * (Rcf)(ub & 1);
    res = RAdd(res, tmp, q);

    ub >>= 1;
    qCtr >>= 1;

    if (qCtr > 0) {
      a = RAdd(a, a, q); // double a mod q for next iteration
    }
  }

  return res; // It is in the correct range
#endif
}

// Quick reversal of a logn-bit dimension
static inline Dim BRev(Dim x, Dim logn) {
  Dim r = x;

  // Reversal by 1 bit, 2 bits, 4 bits, ..., 16 bits
  r = (r & 0x55555555) << 1 | (r & 0xAAAAAAAA) >> 1;
  r = (r & 0x33333333) << 2 | (r & 0xCCCCCCCC) >> 2;
  r = (r & 0x0F0F0F0F) << 4 | (r & 0xF0F0F0F0) >> 4;
  r = (r & 0x00FF00FF) << 8 | (r & 0xFF00FF00) >> 8;
  r = (r & 0x0000FFFF) << 16 | (r & 0xFFFF0000) >> 16;

  // Shift to the size of n
  r = r >> (32 - logn);

  return r;
}

// Some useful constants
constexpr Rcf LWE_CUTOFF = 5;     // 5 standard deviations is the cutoff for SampleZ
constexpr double LWE_EPS = 1e-10; // A double value below LWE_EPS is considered equal to 0

// Parameters for KP- and CP- ABE schemes
constexpr Dim KPABE_L =
    32; // Maximum number of attributes allowed per conjuction for the KP-ABE scheme
constexpr Dim CPABE_L =
    8; // Maximum number of attributes allowed per conjuction for the CP-ABE scheme
constexpr const char *KPABE_GROUPPARAMS = "lweKpabeGroupParams";
constexpr const char *CPABE_GROUPPARAMS = "lweCpabeGroupParams";
