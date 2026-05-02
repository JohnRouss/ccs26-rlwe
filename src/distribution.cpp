#include "distribution.h"

///
/// \file   distribution.cpp
///
/// \brief  Implementation for the one dimensional Gaussian distribution over Z.
///

#include "cpabe.h"
#include "distribution.h"
#include "field_element.h"
#include "kpabe.h"
#include "lwe_base_types.h"
#include "ring.h"
#include "ring_element.h"
#include "trapdoor.h"
#include <cmath>

using namespace std;

/********************************************************************************
 * Implementation of the Distribution class
 ********************************************************************************/

/*!
 * Constructors for the Distribution class.
 *
 */
Distribution::Distribution(RNG *rng) {
  this->m_RNG_ = rng;

  // Low and high should be set by the sub-classes
  this->m_low = 0;
  this->m_high = 0;
}

/*!
 * Destructor for the Distribution class.
 *
 */

Distribution::~Distribution {
}

void Distribution::PrintHistogram(int NumOfSamples, int NumOfBuckets) {
  // Make sure that we have a positive odd number of buckets (so that there is
  // a middle bucket which usually contains 0)
  if (NumOfBuckets < 0) {
    NumOfBuckets = 0;
  }
  NumOfBuckets |= 1;
  std::vector<int> freq(NumOfBuckets, 0);

  // Calculate the bucket size which might exceed the low and high bounds
  Rcf low = this->m_low;
  Rcf high = this->m_high;
  Rcf bucketSz = (high - low + (Rcf)NumOfBuckets) / (Rcf)NumOfBuckets;
  bucketSz |= 1; // Make sure it is odd

  // Center on the mid point and calculate
  // the left edge of the smallest bucket
  Rcf mid = (low + high) / 2;
  Rcf llow = mid;                        // Start from the mid point
  llow -= (NumOfBuckets / 2) * bucketSz; // Go left half buckets
  llow -= (bucketSz / 2);                // Go left half a bucket (for the middle bucket)

  // Sample
  for (int i = 0; i < NumOfSamples; i++) {
    Rcf sample = Sample;
    sample -= llow;
    freq[sample / bucketSz] += 1;
  }

  // Make sure we don't print too long columns.
  // So scale by some amount.
  constexpr int maxC = 30;
  int scale = 1 + NumOfSamples / (NumOfBuckets * maxC);

  // Histogram
  Rcf left = llow;
  Rcf right = llow + bucketSz - 1;
  for (int i = 0; i < NumOfBuckets; i++) {
    freq[i] /= scale;
    cout << "[" << std::setw(5) << left << ", " << std::setw(5) << right << "] ";
    if ((left < mid) && (right > mid)) {
      cout << " ";
    } else {
      cout << "|";
    }
    for (int j = 0; j < freq[i]; j++) {
      cout << "~";
    }
    cout << endl;

    // New edges
    left = right + 1;
    right += bucketSz;
  }
}

ostream &operator<<(ostream &os, Distribution &distr) {
  os << "[low, high) = [" << distr.m_low << ", " << distr.m_high << ")";
  return os;
}

/********************************************************************************
 * Implementation of the UniformDistribution class
 ********************************************************************************/

/*!
 * Constructors for the UniformDistribution class.
 *
 */

UniformDistribution::UniformDistribution(RNG *rng, Rcf q) : Distribution(rng) {
  this->q = q;

  Dim qbits = 0;
  while (q > 0) {
    qbits++;
    q >>= 1;
  }
  this->logq = qbits;

  this->m_high = (this->q / 2) + 1;
  this->m_low = -(this->q / 2);

  this->qBytes = (this->logq + 7) / 8; // Ceiling of number of bytes needed

  int bitsOfMsb = this->logq % 8;         // Number of bits for the most significant byte
  this->msbMask = (~((-1) << bitsOfMsb)); // Zero out the top bits
}

/*!
 * Destructor for the UniformDistribution class.
 *
 */

UniformDistribution::~UniformDistribution {
}

Rcf UniformDistribution::Sample {
  RNG *myRNG = m_RNG_;

  Rcf q = this->q;
  Rcf mask = this->msbMask;
  size_t numOfBytes = this->qBytes;
  Rcf tmp = 0;
  do {
    std::vector<uint8_t> buf(numOfBytes, 0);
    myRNG->getRandomBytes(buf.data, numOfBytes);

    auto it = buf.begin;
    tmp = (Rcf)*it; // Insert the most significant byte
    tmp &= mask;    // Zero out the top bits
    it++;

    while (it != buf.end) {
      tmp = (tmp << 8) | ((Rcf)*it); // shift left 8 bits and add the new byte
      it++;
    }
  } while (tmp >= q); // Reject if bigger than q

  return RRange(tmp, q);
}

ostream &operator<<(ostream &os, UniformDistribution &distr) {
  os << "Uniform distr. in [" << distr.m_low << ", " << distr.m_high << ") w/ (q, logq) = ("
     << distr.q << ", " << distr.logq << ")";
  return os;
}

/********************************************************************************
 * Implementation of the GaussianIntDistribution class
 ********************************************************************************/

/*!
 * Constructor for the GaussianIntDistribution class.
 *
 */

GaussianIntDistribution::GaussianIntDistribution(RNG *rng, Fcf c, Fcf sigma)
    : Distribution(rng), c{c}, sigma{sigma} {

  // Maximum and minimum integer values that the distribution
  // can take. Never exceeding LWE_CUTOFF * sigma distance.
  this->m_high = (Rcf)(std::floor(c + LWE_CUTOFF * sigma));
  this->m_low = (Rcf)(std::ceil(c - LWE_CUTOFF * sigma));

  // Create a seed for the generator
  RNG *myRNG = m_RNG_;
  std::vector<uint8_t> buf(sizeof(uint64_t), 0);
  myRNG->getRandomBytes(buf.data, sizeof(uint64_t));

  uint64_t seed = 0;
  for (uint8_t b : buf) {
    seed = (seed << 8) | b; // shift left 8 bits and add the new byte
  }

  // Initialize the generator and the normal distribution
  this->generator = std::make_unique<mt19937>(seed);
  this->normal = std::make_unique<normal_distribution<Fcf>>(c, sigma);
}

/*!
 * Destructor for the GaussianIntDistribution class.
 *
 */

GaussianIntDistribution::~GaussianIntDistribution {
}

// For sampling we use the C++ standard function from <random>.
// This will provide a Fcf (non integer) number X
// with mean c and standard deviation s.
Rcf GaussianIntDistribution::Sample {
  Rcf high = this->m_high;
  Rcf low = this->m_low;

  Rcf resInt = (Rcf)0;
  do {
    Fcf res = (*this->normal)(*this->generator);
    resInt = (Rcf)round(res);
  } while ((resInt > high) || (resInt < low));

  return resInt;
}

ostream &operator<<(ostream &os, GaussianIntDistribution &distr) {
  os << "GaussianInt distr. in [" << distr.m_low << ", " << distr.m_high << ") w/ (c, sigma) = ("
     << distr.c << ", " << distr.sigma << ")";
  return os;
}

/********************************************************************************
 * Implementation of the GaussianErfDistribution class
 ********************************************************************************/

// For sampling we use the C++ erf function for cmath.
inline Fcf GaussianCdf(Fcf x, Fcf c, Fcf sigmaSqrt2) {
  return 0.5 * (1.0 + std::erf((x - c) / sigmaSqrt2));
}

/*!
 * Constructor for the GaussianErfDistribution class.
 *
 */
GaussianErfDistribution::GaussianErfDistribution(RNG *rng, Fcf c, Fcf sigma)
    : Distribution(rng), sigma{sigma} {

  // Calculate the integer and fractional parts of c
  Fcf doubleIc = std::floor(c);
  this->ic = (Rcf)doubleIc;
  this->fc = c - doubleIc;

  // Approximate maximum and minimum integer values that the distribution
  // can take centered on the fractional part of the center
  this->m_flow = (Rcf)(std::ceil(fc - LWE_CUTOFF * sigma));
  this->m_fhigh = (Rcf)(std::floor(fc + LWE_CUTOFF * sigma));

  // Shift the flow and fhigh
  this->m_low = this->m_flow + this->ic;
  this->m_high = this->m_fhigh + this->ic;

  // Count the number of range bits
  Rcf range = this->m_fhigh - this->m_flow;
  Dim bits = 0;
  while (range > 0) {
    range >>= 1;
    bits++;
  }
  requireCondition(bits <= GaussianErfDistributionERF_MAX_ITERATIONS, statuscode_rand_insufficient);
  this->rangeBits = bits;

  this->sigmaSqrt2 = sigma * std::sqrt(2.0);

  // Add the probabilities of low and high in the table
  // (So that we exclude anything out of this range)
  this->cdfProbs[m_flow] = GaussianCdf((Fcf)m_flow, this->fc, this->sigmaSqrt2);
  this->cdfProbs[m_fhigh] = GaussianCdf((Fcf)m_fhigh, this->fc, this->sigmaSqrt2);

  // Compute the probabilies of the first levels of the binary tree
  // Algorithm: Put pairs of endpoints into a queue
  int precompLevels = GaussianErfDistributionERF_PRECOMPUTED_LVLS;
  std::queue<std::pair<Rcf, Rcf>> endpoints;
  endpoints.push({this->m_flow, this->m_fhigh});
  while (precompLevels > 0) {
    std::pair<Rcf, Rcf> curr;
    do {
      curr = endpoints.front;

      // Find the middle point
      Rcf left = curr.first;
      Rcf right = curr.second;
      Rcf mid = left + right;
      mid >>= 1;

      // Add the cumulative probability in the table
      this->cdfProbs[mid] = GaussianCdf((Fcf)mid, this->fc, this->sigmaSqrt2);

      // Add two new points in the queue
      endpoints.push({left, mid});
      endpoints.push({mid, right});

      endpoints.pop;
    } while (curr.second != this->m_fhigh); // Break if it is the end of the level

    precompLevels--;
  }
}

// Get a random double between [0,1)
Fcf Random0_1(RNG *rng) {
  std::vector<uint8_t> buf(sizeof(Rcf), 0);
  rng->getRandomBytes(buf.data, sizeof(Rcf));

  // Set the sign bit to +
  buf[0] = buf[0] & 0x7F;

  Rcf x = 0;
  for (uint8_t b : buf) {
    x = (x << 8) | b; // shift left 8 bits and add the new byte
  }

  return ((Fcf)x / (Fcf)(RCF_MAX));
}

// Inversion method for sampling
// Without the full precomputed table
Rcf GaussianErfDistribution::Sample {
  // Initial values of low and high indexes
  Rcf low = this->m_flow;
  Rcf high = this->m_fhigh;

  // Pick a random double which will give
  // a value between m_low and m_high
  Fcf lowCdf = this->cdfProbs[low];
  Fcf highCdf = this->cdfProbs[high];
  Fcf U = 0;
  do {
    U = Random0_1(this->m_RNG_);
  } while ((U < lowCdf) || (U > highCdf));

  // Current value of mid
  Rcf mid = 0;

  // Total number of iterations (constant per sigma)
  Dim iterations = this->rangeBits;
  do {
    //
    // The binary search works by calculating the CDF
    // of the mid-point between low and high.
    // If the random probability U is:
    //      - Strictly less than P(mid): change high to mid.
    //      - More or equal to P(mid): change low to mid
    //
    mid = low + high;
    mid = mid >> 1;

    Fcf midProb;
    if (this->rangeBits - iterations < GaussianErfDistributionERF_PRECOMPUTED_LVLS) {
      // For the first GaussianErfDistributionERF_PRECOMPUTED_LVLS we have the
      // probabilities in the precomputed table. Use it.
      midProb = this->cdfProbs[mid];
    } else {
      // Calculate it on the fly
      midProb = GaussianCdf((Fcf)(mid), this->fc, this->sigmaSqrt2);
    }

    // Constant time conditional change
    Rcf mask = (Rcf)(U < midProb);
    mask = mask - 1; // If (U < midProb), then mask = 0, else mask = 0xFFF...F
    high = (mask & high) | (~mask & mid);
    low = (~mask & low) | (mask & mid);

    iterations--;
  } while (iterations > 0);

  return this->ic + low;
}

ostream &operator<<(ostream &os, GaussianErfDistribution &distr) {
  os << "GaussianErf distr. in [" << distr.m_low << ", " << distr.m_high << ") w/ (c, sigma) = ("
     << distr.ic << ", " << distr.sigma << ")";
  return os;
}
