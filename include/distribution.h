///
/// \file   distribution.h
///
/// \brief  Class definition file for different distributions.
///

#include <random>

/// \class    Distribution
/// \brief    Generic container for distributions over Z or Z_q
///
class Distribution : public ObjectBase {
public:
  Distribution(RNG *rng);
  ~Distribution;

  void PrintHistogram(int NumOfSamples, int NumOfBuckets);

  virtual Rcf Sample = 0;

  friend std::ostream &operator<<(std::ostream &os, Distribution &distr);

protected:
  RNG *m_RNG_;
  Rcf m_low;  // Closed
  Rcf m_high; // Open (the sample can go up to high - 1)
};

/// \class    UniformDistribution
/// \brief    Container for the uniform distribution over Z_q
class UniformDistribution : public Distribution {
public:
  UniformDistribution(RNG *rng, Rcf q);
  ~UniformDistribution;

  Rcf Sample;

  friend std::ostream &operator<<(std::ostream &os, UniformDistribution &distr);

protected:
  Rcf q;         // Modulus q
  Dim logq;      // Numbers of bits of q
  size_t qBytes; // Number of bytes of q
  Rcf msbMask;   // Bit-mask to zero out the top bits of the MSB of q
};

/// \class  GaussianIntDistribution
/// \brief  Generic container for the one dimensional Gaussian distribution over Z.
///         It takes as inputs the center of the distribution c,
///         and the standard deviation sigma.
class GaussianIntDistribution : public Distribution {
public:
  GaussianIntDistribution(RNG *rng, Fcf c, Fcf sigma);
  ~GaussianIntDistribution;

  Rcf Sample;

  friend std::ostream &operator<<(std::ostream &os, GaussianIntDistribution &distr);

protected:
  Fcf c;
  Fcf sigma;

  std::unique_ptr<std::mt19937> generator;
  std::unique_ptr<std::normal_distribution<Fcf>> normal;
};

/// \class  GaussianErfDistribution
/// \brief  Container for the one dimensional Gaussian distribution over Z
///         using the C++ erf function.
///         It takes as inputs the center of the distribution c,
///         and the standard deviation sigma.
class GaussianErfDistribution : public Distribution {
public:
  GaussianErfDistribution(RNG *rng, Fcf c, Fcf sigma);
  ~GaussianErfDistribution = default;

  Rcf Sample;

  friend std::ostream &operator<<(std::ostream &os, GaussianErfDistribution &distr);

protected:
  Rcf ic; // Integer part of the center c

  Fcf fc;         // Fractional part of the center c
  Fcf sigma;      // Standard deviation
  Fcf sigmaSqrt2; // S.Dev. multiplied by sqrt(2)

  Rcf m_flow;    // m_low centered on fc
  Rcf m_fhigh;   // m_high centered on fc
  Dim rangeBits; // Number of bits for the range between low and high (to restrict binary search)

  std::map<Rcf, Fcf> cdfProbs; // Precomputed probabilities
};

#define GaussianErfDistributionERF_MAX_ITERATIONS                                                  \
  (50) // Maximum number of binary search iterations (it should not be
       // more than the number of bits of the fraction part of Fcf)
#define GaussianErfDistributionERF_PRECOMPUTED_LVLS                                                \
  (3) // How many levels of the binary tree are precomputed

#ifdef __RCF64__
#define GaussianErfDistribution                                                                    \
  GaussianErfDistribution // Class to be used for all the Gaussians
#else
#define GaussianErfDistribution                                                                    \
  GaussianIntDistribution // Class to be used for all the Gaussians
#endif
