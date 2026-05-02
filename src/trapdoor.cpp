#include "trapdoor.h"

///
/// \file   trapdoor.cpp
///
/// \brief  Implementation for the trapdoor generation and preimage sampling
///         algorithms for Ring LWE
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
 * Implementation of the Trapdoor class
 ********************************************************************************/

/*!
 * Constructors for the Trapdoor class.
 *
 */
Trapdoor::Trapdoor(RNG *rng, std::shared_ptr<Ring> ring) {
  this->m_RNG_ = rng;
  this->ring = ring;
}

/*!
 * Destructor for the Trapdoor class.
 *
 */

Trapdoor::~Trapdoor {
}

/********************************************************************************
 * The following implementations of TrapGen, SampleG, SamplePz, etc.
 * are adapted from:
 *     "Implementation and Evaluation of Improved Gaussian Sampling for Lattice Trapdoors"
 *     by K. Gur, Y. Polyakov, K. Rohloff, G. Ryan, and E. Savas
 * and
 *      "Faster Gaussian Sampling for Trapdoor Lattices with Arbitrary Modulus"
 *      by N. Genise, and D. Micciancio
 ********************************************************************************/

StatusCode Trapdoor::TrapGen(double sigma, std::vector<RingElement> &A,
                             std::vector<RingElement> &T) {
  StatusCode err = STATUS_OK;

  // Get a uniformly random a
  RingElement a(this->ring);
  a.setUniform(this->getRNG, /*timeDomain = */ false);

  // Create a Gaussian distribution centered on 0
  GaussianErfDistribution distr(this->getRNG, 0, sigma);

  // Clear the output vectors
  A.clear;
  T.clear;

  // Create some temporary and constant ring elements
  RingElement two(this->ring);
  two = 2;
  two.ToFrequencyD;

  RingElement pow(this->ring);
  pow = 1;
  pow.ToFrequencyD;

  RingElement sum(this->ring);
  RingElement tmp(this->ring);

  // Add the first two elements to A
  A.push_back(pow);
  A.push_back(a);

  // Main loop for k iterations:
  //  - Sample vector T with 2*k elemenents
  //    The elements r and e will be interleaved in T:
  //    [ r1, e1, r2, e2, ..., rk, ek]
  //  - Calculate the powers of 2
  //  - Store the sum into a
  Rcf k = this->ring->logq;
  for (Rcf i = 0; i < k; i++) {
    // 2^{i}
    sum = pow;

    // get ri
    tmp.setRandom(&distr);
    tmp.ToFrequencyD;
    T.push_back(tmp);

    // -a*ri
    tmp *= a;
    sum -= tmp;

    // get ei
    tmp.setRandom(&distr);
    tmp.ToFrequencyD;
    T.push_back(tmp);

    // -ei
    sum -= tmp;

    A.push_back(sum);

    pow *= two;
  }

  // Error handling omitted in pseudocode artifact.

  return err;
}

//
// Helper functions for SampleG
//

// The following retuns the least significant bit
// of an integer v >= 0 and mutates the value for
// the next bit.
static inline double lsbAndMove(Rcf &v) {
  double res = (double)(v & 1);
  v >>= 1;
  return res;
}

vector<Rcf> SampleD(RNG *rng, double sigma, vector<double> &c, vector<double> &d) {
  requireCondition(c.size == d.size, statuscode_invalid_length);

  // Get the adjusted c vector
  GaussianErfDistribution d1(rng, -c.back / d.back, sigma / d.back);
  Rcf zk_1 = d1.Sample;

  vector<double> cPrime;
  for (int i = 0; i < c.size; i++) {

    // Note: There is probably a sign typo in the following line
    // in the "Improved ..." paper (See SampleD in  "Faster ...").
    double tmpD = c[i] + (double)zk_1 * d[i];
    cPrime.push_back(tmpD);
  }

  // Sample the other elements
  vector<Rcf> z;
  for (int i = 0; i < c.size - 1; i++) {
    GaussianErfDistribution d2(rng, -cPrime[i], sigma);
    Rcf tmpZ = d2.Sample;
    z.push_back(tmpZ);
  }

  // Add the last element
  z.push_back(zk_1);

  return z;
}

vector<Rcf> Perturb(RNG *rng, Dim k, double sigma, vector<double> &l, vector<double> &h) {
  requireCondition(l.size == h.size, statuscode_invalid_length);

  double beta = 0;
  vector<Rcf> z;
  for (Dim i = 0; i < k; i++) {
    GaussianErfDistribution distr(rng, beta / l[i], sigma / l[i]);
    Rcf tmpZ = distr.Sample;
    z.push_back(tmpZ);

    beta = -tmpZ * h[i];
  }

  vector<Rcf> p;
  Rcf tmpP = 5 * z[0] + 2 * z[1];
  p.push_back(tmpP);
  for (Dim i = 1; i < k - 1; i++) {
    tmpP = 2 * (z[i - 1] + 2 * z[i] + z[i + 1]);
    p.push_back(tmpP);
  }
  tmpP = 2 * (z[k - 2] + 2 * z[k - 1]);
  p.push_back(tmpP);

  return p;
}

vector<RingElement> Trapdoor::SampleG(double s, RingElement u) {
  RNG *rng = this->getRNG;
  Rcf q = this->ring->q;
  Dim n = this->ring->n;
  Dim logq = this->ring->logq;

  double sigma = s / 3;
  double k = (double)logq;

  Rcf qtmp = q;

  // Possible Optimization: Vectors l, h, d can be precomputed for fixed q.
  // Consider this to marginally improve the performance of KeyGen, where
  // SampleG is called once.
  vector<double> l, h, d;

  l.push_back(sqrt(2 * (1 + 1 / k) + 1));
  h.push_back(0);

  double dtmp = (lsbAndMove(qtmp)) / 2;
  d.push_back(dtmp);

  for (Dim i = 1; i < logq; i++) {
    l.push_back(sqrt(2 * (1 + 1 / (k - (double)i))));
    h.push_back(sqrt(2 * (1 - 1 / (k - (double)i + 1))));

    dtmp = (dtmp + lsbAndMove(qtmp)) / 2;
    d.push_back(dtmp);
  }

  vector<Rcf> t; // This vector will hold all the k x n elements
  for (Dim i = 0; i < n; i++) {
    Rcf v = u[i];
    v += q * (Rcf)(v < 0); // To get the correct bits of v move it to the [0,q-1] range
    vector<Rcf> p = Perturb(rng, logq, sigma, l, h);

    vector<double> c;
    double tmpC = (lsbAndMove(v) - (double)p[0]) / 2;
    c.push_back(tmpC);

    for (Dim j = 1; j < logq; j++) {
      v >>= 1;
      tmpC = (tmpC + (lsbAndMove(v) - (double)p[j])) / 2;
      c.push_back(tmpC);
    }

    vector<Rcf> z = SampleD(rng, sigma, c, d);

    v = u[i];
    v += q * (Rcf)(v < 0); // To get the correct bits of v move it to the [0,q-1] range
    qtmp = q;
    Rcf tmpT = 2 * z[0] + lsbAndMove(qtmp) * z.back + lsbAndMove(v);
    t.push_back(tmpT);

    for (Dim j = 1; j < logq - 1; j++) {
      tmpT = 2 * z[j] + lsbAndMove(qtmp) * z.back + lsbAndMove(v) - z[j - 1];
      t.push_back(tmpT);
    }

    tmpT = lsbAndMove(qtmp) * z.back + lsbAndMove(v) - z[logq - 2];
    t.push_back(tmpT);
  }

  // Convert to ring elements
  vector<RingElement> res;
  for (Dim i = 0; i < logq; i++) {
    vector<Rcf> zTmp;
    RingElement tmp(this->ring);
    for (Rcf j = 0; j < n; j++) {
      // Take the elements every logq indexes
      Rcf tmpT = t[j * logq + i] % q; // This can give a negative number with abs < q
      zTmp.push_back(RRange(tmpT, q));
    }
    tmp.setValues(zTmp, /*timeD = */ true);
    tmp.ToFrequencyD;
    res.push_back(tmp);
  }

  return res;
}

//
// Helper functions for SamplePz
//

vector<Rcf> SampleFz(RNG *rng, FieldElement &f, FieldElement &c);

vector<Rcf> Sample2z(RNG *rng, FieldElement &a, FieldElement &b, FieldElement &d, FieldElement &c0,
                     FieldElement &c1) {

  vector<Rcf> q1 = SampleFz(rng, d, c1);

  FieldElement q1Hat = ConvertIntegersToField(q1);

  // Field math in frequency domain
  FieldSpectrum q1HatS = Fft(q1Hat);
  FieldSpectrum c0S = Fft(c0);
  FieldSpectrum c1S = Fft(c1);
  FieldSpectrum aS = Fft(a);
  FieldSpectrum bS = Fft(b);
  FieldSpectrum dS = Fft(d);

  // Main operations
  FieldSpectrum bDivD = bS * ~dS; // b / d
  c0S += bDivD * (q1HatS - c1S);  // c0 = c0 + b/d * (q1 - c1)
  aS -= bDivD * !bS;              // a = a - b/d * b^T

  // Convert to time domain
  FieldElement c0New = Ifft(c0S);
  FieldElement fNew = Ifft(aS);

  vector<Rcf> q0 = SampleFz(rng, fNew, c0New);

  // Concatenate the two vectors (q0, q1)
  vector<Rcf> qRes;
  qRes.reserve(q0.size + q1.size);
  qRes.insert(qRes.end, q0.begin, q0.end);
  qRes.insert(qRes.end, q1.begin, q1.end);
  return qRes;
}

vector<Rcf> SampleFz(RNG *rng, FieldElement &f, FieldElement &c) {
  vector<Rcf> res;
  if (f.size == 1) {
    requireCondition(c.size == 1, statuscode_invalid_length);
    requireCondition(f[0] > 0, statuscode_invalid_field_element_sign);

    GaussianErfDistribution distr(rng, c[0], sqrt(f[0]));

    // Get a sample from the one dimensional Gaussian distribution.
    // Note: There is no need to mod the result into the [-q/2, q/2]
    // range because its absolute value should be smaller than q/4.
    Rcf tmp = distr.Sample;

    res.push_back(tmp);
  } else {
    Dim n2 = ((Dim)f.size) >> 1;

    FieldElement f0(n2);
    FieldElement f1(n2);
    ExtractEvenAndOdd(f, f0, f1);

    FieldElement c0(n2);
    FieldElement c1(n2);
    ExtractEvenAndOdd(c, c0, c1);

    vector<Rcf> qRes = Sample2z(rng, f0, f1, f0, c0, c1);

    // Interleave q0 and q1 from Sample2z
    requireCondition(qRes.size == f.size, statuscode_invalid_length);
    for (Dim i = 0; i < n2; i++) {
      res.push_back(qRes[i]);
      res.push_back(qRes[i + n2]);
    }
  }
  return res;
}

vector<RingElement> Trapdoor::SamplePz(double s, double alpha, vector<RingElement> &T) {
  RNG *rng = this->getRNG;

  double sSqr = s * s;
  double alphaSqr = alpha * alpha;

  // Calculate z = 1 / (1/alpha^2 - 1/s^2)
  double z = 1.0 / (1.0 / alphaSqr - 1.0 / sSqr);

  Dim k = this->ring->logq;
  Dim n = this->ring->n;

  // Initialize a, b, d to zero (field) elements
  FieldElement a(n), b(n), d(n);

  // Construct the sums for a, b, d in the same loop
  requireCondition(T.size == this->sizeOfT, statuscode_invalid_length);
  for (Dim i = 0; i < k; i++) {
    RingElement riSq = !T[2 * i] * T[2 * i];         // ri^t * ri
    RingElement riei = !T[2 * i] * T[2 * i + 1];     // ri^t * ei
    RingElement eiSq = !T[2 * i + 1] * T[2 * i + 1]; // ei^t * ri

    // Convert them to time domain
    riSq.ToTimeD;
    riei.ToTimeD;
    eiSq.ToTimeD;

    // Convert them to field elements
    FieldElement riSqF = ConvertIntegersToField(riSq);
    FieldElement rieiF = ConvertIntegersToField(riei);
    FieldElement eiSqF = ConvertIntegersToField(eiSq);

    // Add them into a, b, d
    a += riSqF;
    b += rieiF;
    d += eiSqF;
  }

  // Scale the sums by z
  a *= z;
  b *= z;
  d *= z;

  // Calculate b = - z Sum(ri^t * ei), a = s^2 - z Sum(ri^t * ri),
  // d = s^2 - z Sum(ei^t * ei)
  FieldElement sSrqF(n); // sSrqF = 0
  b = sSrqF - b;
  sSrqF = sSqr; // sSrqF = s^2
  a = sSrqF - a;
  d = sSrqF - d;

  // Calculate the k q_i elements
  // using the SampleZ distribution with
  // standard deviation sqrt(s^2 - alpha^2)
  // In the same loop construct c0 and c1 as ring elements where
  // c0R = Sum(r_i * q_i) and c1R = Sum(e_i * q_i)
  double stdev = std::sqrt(sSqr - alphaSqr);
  GaussianErfDistribution distr(rng, 0, stdev);

  RingElement c0R(this->ring, /*timeDomain =*/false); // c0R = 0
  RingElement c1R(this->ring, /*timeDomain =*/false); // c1R = 0

  vector<RingElement> qElems;
  RingElement tmp(this->ring);
  RingElement tmpProd(this->ring);
  for (Dim i = 0; i < k; i++) {
    tmp.setRandom(&distr, /*timeDomain =*/true); // It should be Gaussian in the time domain
    tmp.ToFrequencyD;

    qElems.push_back(tmp);

    // r_i * q_i
    tmpProd = T[2 * i] * tmp;
    c0R += tmpProd;

    // e_i * q_i
    tmpProd = T[2 * i + 1] * tmp;
    c1R += tmpProd;
  }

  // Convert to field elements and scale
  // Note: There is probably a sign typo in the following line
  // in the "Improved ..." paper (See SamplePz in  "Faster ...").
  double factor = -alphaSqr / (sSqr - alphaSqr);
  c0R.ToTimeD;
  c1R.ToTimeD;

  FieldElement c0 = ConvertIntegersToField(c0R);
  FieldElement c1 = ConvertIntegersToField(c1R);

  c0 *= factor;
  c1 *= factor;

  vector<Rcf> p = Sample2z(rng, a, b, d, c0, c1); // p is of size 2*n
  requireCondition(p.size == 2 * n, statuscode_invalid_length);

  vector<RingElement> res;

  // Construct and insert the first two elements
  RingElement p0(this->ring);
  vector<Rcf> p0V(p.begin, p.begin + n);
  p0.setValues(p0V);
  p0.ToFrequencyD;
  res.push_back(p0);

  RingElement p1(this->ring);
  vector<Rcf> p1V(p.begin + n, p.end);
  p1.setValues(p1V);
  p1.ToFrequencyD;

  res.push_back(p1);

  // After that insert all the q elements
  for (RingElement rElem : qElems) {
    res.push_back(rElem);
  }

  return res;
}

// Helper functions to calculate the s and alpha parameters
inline double sParameter(Dim k, Dim n, double sigma) {
  // C is selected empirically larger than the value commonly cited
  // in the literature (~1.8) so that the covariance is not
  // degenerate for the parameter regimes used here.
  constexpr double C = 6.5;
  constexpr double t = 4.7;

  return C * sigma * sigma * (sqrt((double)n * (double)k) + sqrt((double)n * 2.0) + t);
}

inline double alphaParameter(double sigma) {
  return 2.0 * sigma;
}

StatusCode Trapdoor::GaussianPreimageSampling(std::vector<RingElement> &A,
                                              std::vector<RingElement> &T, RingElement &u,
                                              double sigma, std::vector<RingElement> &X) {
  StatusCode err = STATUS_OK;

  X.clear;
  Dim k = this->ring->logq;
  Dim n = this->ring->n;

  double s = sParameter(k, n, sigma);
  double alpha = alphaParameter(sigma);
  vector<RingElement> p = this->SamplePz(s, alpha, T);

  RingElement sum(ring);
  sum = 0;
  sum.ToFrequencyD;

  requireCondition(p.size == this->sizeOfA, statuscode_invalid_length);
  requireCondition(A.size == this->sizeOfA, statuscode_invalid_length);
  for (Dim i = 0; i < p.size; i++) {
    sum += (p[i] * A[i]);
  }

  // u should be in the time domain
  RingElement uT = u;
  uT.ToTimeD;

  sum.ToTimeD;
  RingElement uPrime(ring);
  uPrime = uT - sum;
  vector<RingElement> z = this->SampleG(sigma, uPrime);

  // First calculate X0 and X1
  RingElement X0 = p[0];
  RingElement X1 = p[1];
  for (Dim i = 0; i < k; i++) {
    RingElement prod1 = z[i] * T[2 * i];     // zi * ri
    RingElement prod0 = z[i] * T[2 * i + 1]; // zi * ei

    X0 += prod0;
    X1 += prod1;
  }
  X.push_back(X0);
  X.push_back(X1);

  for (Dim i = 0; i < k; i++) {
    RingElement tmp = p[i + 2] + z[i];
    X.push_back(tmp);
  }

  // Error handling omitted in pseudocode artifact.

  return err;
}

Dim Trapdoor::sizeOfA {
  return 2 + this->ring->logq; // k + 2
}

Dim Trapdoor::sizeOfT {
  return 2 * this->ring->logq; // 2 * k
}

double Trapdoor::sampleStdev(double sigma) {
  Dim k = this->ring->logq;
  Dim n = this->ring->n;

  double sSqr = sParameter(k, n, sigma);
  sSqr *= sSqr;

  double alphaSqr = alphaParameter(sigma);
  alphaSqr *= alphaSqr;

  return std::sqrt(sSqr - alphaSqr);
}