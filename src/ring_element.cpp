#include "ring_element.h"

///
///	\file   ring_element.cpp
///
///	\brief  Class implementation for ring elements.
///

#include <algorithm>
#include <ostream>
#include <vector>

#include "cpabe.h"
#include "distribution.h"
#include "field_element.h"
#include "kpabe.h"
#include "lwe_base_types.h"
#include "ring.h"
#include "ring_element.h"
#include "trapdoor.h"

using namespace std;

/********************************************************************************
 * Implementation of the RingElement class
 ********************************************************************************/
/*!
 * Constructor for the RingElement class.
 *
 */
RingElement::RingElement(std::shared_ptr<Ring> ring, bool timeDomain) {
  if (ring == nullptr) {
    throw statuscode_invalid_group_params;
  }
  this->assign(ring->n, 0);
  this->ring = ring;
  this->timeD = timeDomain;
}

std::shared_ptr<Ring> RingElement::getRing const {
  return this->ring;
}

bool RingElement::isTimeDomain const {
  return this->timeD;
}

void RingElement::setValues(std::vector<Rcf> &input, bool timeDomain) {
  Dim n = this->ring->n;
  size_t count = std::min<size_t>(n, input.size);

  Rcf q = this->ring->q;
  for (Dim i = 0; i < count; i++) {
    (*this)[i] = RRange(input[i], q);
  }

  // Fill the rest with zeros (if needed)
  for (Dim i = count; i < n; i++) {
    (*this)[i] = 0;
  }

  this->timeD = timeDomain;
}

void RingElement::setConstantTerm(Rcf c, bool timeDomain) {
  Dim n = this->ring->n;
  Rcf q = this->ring->q;

  (*this)[0] = RRange(c, q);
  for (Dim i = 1; i < n; i++) {
    (*this)[i] = 0;
  }

  this->timeD = timeDomain;
}

void RingElement::setRandom(Distribution *distr, bool timeDomain) {
  Dim n = this->ring->n;
  Rcf q = this->ring->q;

  for (Dim i = 0; i < n; i++) {
    (*this)[i] = RRange(distr->Sample, q);
  }

  this->timeD = timeDomain;
}

void RingElement::setUniform(RNG *rng, bool timeDomain) {
  UniformDistribution distr(rng, this->ring->q);
  this->setRandom(&distr, timeDomain);
}

void RingElement::encodeOneZero(std::vector<uint8_t> &msg) {
  Dim n = this->ring->n;

  Rcf q = this->ring->q;
  Rcf q2 = q >> 1; // Floor of q / 2

  this->clear;

  // Get each byte from msg
  // Make sure we don't add more than n components
  vector<uint8_t>::iterator itMsg = msg.begin;
  Dim componentsAdded = 0;
  uint8_t byte = 0;
  while (componentsAdded < n) {
    if (componentsAdded % 8 == 0) {
      if (itMsg != msg.end) {
        byte = *itMsg;
        itMsg++;
      } else {
        byte = 0;
      }
    }

    Rcf tmp = q2 * (Rcf)(byte & 0x1); // Either 0 or q2
    this->push_back(tmp);

    componentsAdded++;
    byte >>= 1;
  }

  this->timeD = true;
}

std::vector<uint8_t> RingElement::decodeOneZero(uint32_t keyByteLen) {
  // Convert to time domain if not there
  if (!this->timeD) {
    this->ToTimeD;
  }

  std::vector<uint8_t> res;

  Rcf q = this->ring->q;
  Rcf bound = q / 4; // Floor of q / 4

  uint32_t bytesDecoded = 0;
  uint8_t tmpByte = 0;

  // ctr indicates how many bits are inserted and
  // how many bits to shift up the current bit
  int ctr = 0;
  for (auto it = this->begin; it != this->end;) {
    uint8_t bitSet = (uint8_t)(*it > bound) | (uint8_t)(*it < -bound);
    bitSet <<= ctr;
    tmpByte |= bitSet;

    ctr++;
    it++; // Increase the iterator here to check if we reached the end

    if ((ctr % 8 == 0) || (it == this->end)) {
      res.push_back(tmpByte);
      tmpByte = 0;
      ctr = 0;
      if (++bytesDecoded >= keyByteLen) {
        break;
      }
    }
  }

  return res;
}

/********************************************************************************
 * The following implementations of NaiveNtt, NaiveIntt, FastNtt, and FastIntt
 * are adapted from:
 *     "A Complete Beginner Guide to the Number Theoretic Transform (NTT)"
 *     by A. Satriawan, R. Mareta, and H. Lee
 ********************************************************************************/

void RingElement::NaiveNtt {
  if (!this->timeD) {
    throw statuscode_invalid_ring_element_domain;
    return;
  }

  if (this->ring->psi == 0) {
    throw statuscode_class_not_initialized;
    return;
  }

  Dim n = this->ring->n;
  Rcf psi = this->ring->psi;
  Rcf q = this->ring->q;

  // Naive algorithm to calculate NTT
  // from the formula a_j = Sum_i (psi^{2ij+i} * a_i).
  // Each term a_j has a multiplier psi^{2ij+i}. This means
  // that the first row multipliers increase by (x psi) each
  // step, the 2nd row by (x 3 psi), the 3rd row by (x 5psi)
  // etc.
  vector<Rcf> vec;
  Rcf psiSq = RMul(psi, psi, q);
  Rcf psiInc = psi; // The first row increment is psi
  for (Dim j = 0; j < n; j++) {
    Rcf res = 0;
    Rcf mult = 1; // The first multiplier is always 1

    for (Rcf ai : *this) {
      Rcf prod = RMul(mult, ai, q);
      res = RAdd(res, prod, q);
      mult = RMul(mult, psiInc, q);
    }

    vec.push_back(res);

    // Increase the psiInc for the next row
    psiInc = RMul(psiInc, psiSq, q);
  }

  // Copy the new vector with bit reversed order
  // and specify the frequency domain
  Dim logn = this->ring->logn;
  for (Dim i = 0; i < n; i++) {
    (*this)[i] = vec[BRev(i, logn)];
  }
  this->timeD = false;
}

void RingElement::NaiveIntt {
  if (this->timeD) {
    throw statuscode_invalid_ring_element_domain;
    return;
  }

  if (this->ring->psiInv == 0) {
    throw statuscode_class_not_initialized;
    return;
  }

  Dim n = this->ring->n;
  Rcf nInv = this->ring->nInv;
  Rcf psiInv = this->ring->psiInv;
  Rcf q = this->ring->q;

  // Naive algorithm to calculate INTT
  // from the formula a_i = n^-1 Sum_j (psi^-{2ij+j} * a_j).
  // Each term a_j has a multiplier psi^{2ij+j}. This means
  // that the first row multipliers start at psi^-0 and
  // increase by (x 1) each step, the 2nd row start at psi^-1
  // and increaseby (x 2 psi), the 3rd row start at psi^-2 by (x 4psi)
  // etc.
  vector<Rcf> vec;
  Dim logn = this->ring->logn;
  Rcf psiSq = RMul(psiInv, psiInv, q);
  Rcf psiInc = 1;   // The first row increment is 1
  Rcf psiFirst = 1; // The first row first term is 1
  for (Dim i = 0; i < n; i++) {
    Rcf res = 0;
    Rcf mult = psiFirst;

    // Note: *this is stored in bit-reversed order
    for (Dim j = 0; j < n; j++) {
      Rcf aj = (*this)[BRev(j, logn)];
      Rcf prod = RMul(mult, aj, q);
      res = RAdd(res, prod, q);
      mult = RMul(mult, psiInc, q);
    }

    // Divide by n
    res = RMul(res, nInv, q);

    vec.push_back((Rcf)res);

    // Increase the psiInc for the next row
    psiInc = RMul(psiInc, psiSq, q);

    // And the psiFirst
    psiFirst = RMul(psiFirst, psiInv, q);
  }

  // Copy the new vector
  // and specify the time domain
  this->clear;
  this->insert(this->end, vec.begin, vec.end);
  this->timeD = true;
}

// Helper function to implement the CT butterfly.
// It takes in the indices of the terms to add
// and the twiddle factor.
// It replaces them in the vector with
// v[i] + twf * v[j] and v[i] - twf * v[j] (Mod q)
void CtButterfly(vector<Rcf> &v, Dim i, Dim j, Rcf twf, Rcf q) {
  Rcf secondTerm = RMul(v[j], twf, q);
  Rcf b1 = RAdd(v[i], secondTerm, q);
  Rcf b2 = RSub(v[i], secondTerm, q);

  v[i] = b1;
  v[j] = b2;
}

void RingElement::FastNtt {
  Rcf q = this->ring->q;
  Dim n = this->ring->n;

  // Assumption: n is a power of 2

  // Check the time domain
  if (!this->timeD) {
    throw statuscode_invalid_ring_element_domain;
    return;
  }

  // jump is how far the other branch of the butterfly is and
  // at the same time the amount of indexes to skip when we
  // switch to a sub-FNTT.
  // For example, for n = 8, the second stage of butterflies
  // does the indices (0, 2), (1, 3), (4, 6), (5, 7).
  // In this case, jump = 2, and the skip
  // happens between (1, 3) and (4, 6).
  Dim nhalf = n >> 1;
  Dim jump = nhalf;

  // The jump also specifies how many twiddle factors we use
  // in bit-reversed order for each stage. (1 for the first
  // stage, 2 for the second, 4 for the third, etc.)
  // We start with the second twiddle factor on index 1, which
  // on the bit-reversed order array is psi^{n/2}
  Dim twfIdx = 1;

  // Go for log(n) stages by shifting jump 1 right each time
  while (jump > 0) {
    Dim upper = 0; // The starting upper index is always 0

    // Every stage of the CT algorithm does n/2 butterflies
    for (Dim i = 0; i < nhalf;) {
      Dim lower = upper + jump;
      Rcf twf = this->ring->psiVecR[twfIdx];

      CtButterfly(*this, upper, lower, twf, q);

      // Increase upper index by 1
      upper++;

      // Check if the next iteration (i+1) needs to jump to a different sub FNTT
      // which moves the upper index and the twiddle factor
      // Note: this is not triggered in the first stage (jump = n/2)
      if ((++i) % jump == 0) {
        upper += jump;
        twfIdx++;
      }
    }

    jump >>= 1;
  }

  // Specify the frequency domain
  // Note: The output is in bit-reversed order!!!
  this->timeD = false;
}

// Helper function to implement the GS butterfly.
// It takes in the indices of the terms to add
// and the twiddle factor.
// It replaces them in the vector with
// v[i] + v[j] and (v[i] - v[j]) * twf (Mod q)
void GsButterfly(vector<Rcf> &v, Dim i, Dim j, Rcf twf, Rcf q) {
  Rcf b1 = RAdd(v[i], v[j], q);
  Rcf b2 = RSub(v[i], v[j], q);
  b2 = RMul(b2, twf, q);

  v[i] = b1;
  v[j] = b2;
}

void RingElement::FastIntt {
  Rcf q = this->ring->q;
  Dim n = this->ring->n;

  // Assumption: n is a power of 2

  // Check for bit reversed ordering
  if (this->timeD) {
    throw statuscode_invalid_ring_element_domain;
    return;
  }

  // jump is how far the other branch of the butterfly is and
  // at the same time the amount of indexes to skip when we
  // switch to a sub-FNTT.
  // It is the reverse of the Fast NTT
  Dim nhalf = n >> 1;
  Dim jump = 1;

  // The jump also specifies how many twiddle factors we use
  // in bit-reversed order for each stage. (n/2 for the first
  // stage, n/4 for the second, n/8 for the third, etc.)
  // We start with the last twiddle factor on index 1, which
  // on the bit-reversed order array is psi^{-(n-1)}
  Dim twfIdx = n - 1;

  // Go for log(n) stages by shifting jump 1 left each time
  while (jump < n) {
    Dim lower = n - 1; // The starting lower index is always n - 1

    // Every stage of the CT algorithm does n/2 butterflies
    for (Dim i = 0; i < nhalf;) {
      Dim upper = lower - jump;
      Rcf twf = this->ring->psiInvVecR[twfIdx];

      GsButterfly(*this, upper, lower, twf, q);

      // Decrease lower index by 1
      lower--;

      // Check if the next iteration (i+1) needs to jump to a different sub FNTT
      // which moves the lower index and the twiddle factor
      // Note: this is not triggered in the last stage (jump = n/2)
      if ((++i) % jump == 0) {
        lower -= jump;
        twfIdx--;
      }
    }

    jump <<= 1;
  }

  // Divide the final result by n
  Rcf nInv = this->ring->nInv;
  for (Dim i = 0; i < n; i++) {
    (*this)[i] = RMul((*this)[i], nInv, q);
  }

  // Specify the time domain
  this->timeD = true;
}

void RingElement::ToTimeD {
  if (!this->timeD) {
    this->FastIntt;
  }
}

void RingElement::ToFrequencyD {
  if (this->timeD) {
    this->FastNtt;
  }
}

RingElement &RingElement::operator=(const Rcf &rhs) {
  if (this->ring == nullptr) {
    throw statuscode_invalid_group_params;
  }

  this->setConstantTerm(rhs, true); // Always time domain and normal order

  return *this;
}

RingElement &RingElement::operator+=(const RingElement &rhs) {
  if ((this->ring->q != rhs.ring->q) || (this->ring->n != rhs.ring->n) ||
      (this->timeD != rhs.timeD)) {
    throw statuscode_invalid_ring_element_domain;
  }

  Rcf q = this->ring->q;
  for (Dim i = 0; i < this->ring->n; i++) {
    (*this)[i] = RAdd((*this)[i], rhs[i], q);
  }

  return *this;
}

RingElement &RingElement::operator-=(const RingElement &rhs) {
  if ((this->ring->q != rhs.ring->q) || (this->ring->n != rhs.ring->n) ||
      (this->timeD != rhs.timeD)) {
    throw statuscode_invalid_ring_element_domain;
  }

  Rcf q = this->ring->q;
  for (Dim i = 0; i < this->ring->n; i++) {
    (*this)[i] = RSub((*this)[i], rhs[i], q);
  }

  return *this;
}

RingElement componentWiseMul(RingElement &lhs, const RingElement &rhs) {
  if ((lhs.getRing->q != rhs.getRing->q) || (lhs.getRing->n != rhs.getRing->n) ||
      (lhs.isTimeDomain != rhs.isTimeDomain)) {
    throw statuscode_invalid_ring_element_domain;
  }

  Rcf q = lhs.getRing->q;
  Dim n = lhs.getRing->n;
  RingElement z(lhs.getRing, lhs.isTimeDomain);

  for (Dim i = 0; i < n; i++) {
    z[i] = RMul(lhs[i], rhs[i], q);
  }

  return z;
}

RingElement negacyclicNaiveMul(RingElement &lhs, const RingElement &rhs) {
  if ((lhs.getRing->q != rhs.getRing->q) || (lhs.getRing->n != rhs.getRing->n) ||
      (lhs.isTimeDomain != rhs.isTimeDomain)) {
    throw statuscode_invalid_ring_element_domain;
  }

  Rcf q = lhs.getRing->q;
  Dim n = lhs.getRing->n;
  RingElement z(lhs.getRing, lhs.isTimeDomain);

  // First make a reverse copy of the rhs vector
  vector<Rcf> revRhs(rhs.rbegin, rhs.rend);

  // Then go through both vectors multiplying elements and adding (or subtracting)
  // monomials according to the negative wrapped convolution.
  for (Dim k = 0; k < n; k++) {
    Rcf acc = 0;

    // First put the last element of reverse rhs in front
    std::rotate(revRhs.begin, revRhs.end - 1, revRhs.end);

    for (Dim i = 0; i < n; i++) {
      Rcf monomial = RMul(lhs[i], revRhs[i], q);

      if (i > k) {
        // Subtract mod q
        acc = RSub(acc, monomial, q);
      } else {
        // Add mod q
        acc = RAdd(acc, monomial, q);
      }
    }

    z[k] = acc;
  }

  return z;
}

RingElement &RingElement::operator*=(const RingElement &rhs) {
  if (this->timeD != rhs.timeD) {
    throw statuscode_invalid_ring_element_domain;
  } else if (this->timeD) {
    *this = negacyclicNaiveMul(*this, rhs); // O(n^2)
  } else {
    *this = componentWiseMul(*this, rhs); // O(n)
  }
  return *this;
}

RingElement &RingElement::operator*=(const Rcf &scaleR) {
  const Dim n = this->ring->n;
  const Rcf q = this->ring->q;

  for (Dim i = 0; i < n; i++) {
    (*this)[i] = RMul((*this)[i], scaleR, q);
  }

  return *this;
}

RingElement RingElement::operator!{
  RingElement z(this->ring, this->timeD);
  Dim n = this->ring->n;

  if (this->timeD) {
    // Transposition in the time domain (definition)
    z[0] = this->at(0);
    for (Dim i = 1; i < n; i++) {
      z[i] = -(this->at(n - i));
    }
  } else {
    // Transposition in the frequency domain
    // with bit-reversed order.
    // Note: It can be proved that:
    //  - Reflection of the coefficients in the
    //    frequency domain gives the frequency domain
    //    coefficients of the transposed element.
    //  - The reflection permutation commutes with
    //    the bit-reverse permutation.
    for (Dim i = 0; i < n; i++) {
      z[i] = (*this)[n - 1 - i];
    }
  }

  return z;
}

ostream &operator<<(ostream &os, RingElement &elem) {
  bool timeDomain = elem.timeD;
  Dim logn = elem.ring->logn;

  size_t width = decimalLength(elem.ring->q);

  os << (timeDomain ? "[ " : "< ");
  for (Dim i = 0; i < elem.ring->n; i++) {
    Dim idx = timeDomain ? i : BRev(i, logn);
    os << std::setw(width + 1) << std::setfill(' ') << elem[idx] << " ";
  }
  os << (timeDomain ? "]" : ">");
  return os;
}

bool operator==(const RingElement &lhs, const RingElement &rhs) {
  if ((lhs.ring->n != rhs.ring->n) || (lhs.ring->q != rhs.ring->q))
    return false;

  Dim logn = lhs.ring->logn;
  int rc = 0;
  for (size_t i = 0; i < lhs.size; i++) {
    Dim lIdx = lhs.timeD ? i : BRev(i, logn);
    Dim rIdx = rhs.timeD ? i : BRev(i, logn);
    rc |= (lhs[lIdx] ^ rhs[rIdx]);
  }

  /* 0 => lhs == rhs, > 0 => lhs != rhs */
  return (rc == 0) ? true : false;
}

bool RingElement::isEqual(ObjectBase *z) const {
  RingElement *z1 = dynamic_cast<RingElement *>(z);
  if (z1 != NULL) {
    return *z1 == *this;
  }
  return false;
}

void RingElement::serialize(std::vector<uint8_t> &result) const {
  result.zeroize;

  // insert the type
  result.insertFirstByte(ELEMENT_RING_ELEMENT);

  // convert to time domain
  RingElement t(*this);
  t.ToTimeD;

  // pack each coefficient
  int byteLen = sizeof(Rcf);
  uint8_t b[sizeof(Rcf)] = {0};
  for (Dim i = 0; i < t.ring->n; i++) {
    Rcf curr = t[i];

    for (int j = 0; j < byteLen; j++) {
      b[j] = (curr & 0xFF); // record last byte
      curr >>= 8;           // shift right by 8 bits
    }

    result.appendArray(b, byteLen);
  }
}

Rcf RingElement::InfNorm {
  if (!this->timeD) {
    throw statuscode_invalid_ring_element_domain; // Only norm in the time domain
  }

  Rcf max = 0;
  for (Dim i = 0; i < this->ring->n; i++) {
    Rcf c = (*this)[i];

    // Take the abs
    Rcf neg = (Rcf)(c < 0);      // neg = 1 if c < 0; neg = 0 otherwise
    neg = neg - 1;               // neg = 0 if c < 0; neg = 0xFF..F otherwise
    c = (c & neg) | (-c & ~neg); // Take the abs(c)

    Rcf big = (Rcf)(c > max);
    big = big - 1;
    max = (max & big) | (c & ~big);
  }

  return max;
}

std::vector<RingElement> hashToRingElements(std::shared_ptr<Ring> ring, std::vector<uint8_t> &key,
                                            const std::string &input, Dim m) {

  // Each call to SHA-256 with prefix the key and suffix an index will generate
  // 256 / (sizeof(Rcf) * 8) ring coefficients. We iterate until we get all the
  // m * n required elements.
  int numOfRcfPerHash = (32 + sizeof(Rcf) - 1) / sizeof(Rcf); // Ceiling of 256 / (sizeof(Rcf) * 8)
  Dim n = (int)ring->n;
  Rcf q = ring->q;
  Dim totalNumberOfRcf = (int)m * (int)ring->n;

  vector<Rcf> ringCoefficients;
  Dim idx = 0; // Index to hash into each SHA256 chunk
  std::vector<uint8_t> tmpInput;
  Rcf tmpRcf;
  while (totalNumberOfRcf > 0) {
    tmpInput = key;
    tmpInput += input;
    tmpInput += std::to_string(idx);

    string digest;
    string str = tmpInput.toString;
    sha256(digest, str);

    uint8_t *xstr = (uint8_t *)digest.c_str;

    tmpRcf = 0;
    for (int i = 0; i < (sizeof(Rcf) * numOfRcfPerHash); i++) {
      tmpRcf <<= 8;
      tmpRcf |= (*xstr);

      xstr++;

      if ((i + 1) % sizeof(Rcf) == 0) {
        tmpRcf = RRange(tmpRcf, q);
        ringCoefficients.push_back(tmpRcf);
        tmpRcf = 0;
      }
    }

    totalNumberOfRcf -= numOfRcfPerHash;
    idx++;
  }

  requireCondition(ringCoefficients.size >= (m * n), statuscode_invalid_input);

  // Pack them all into ring elements
  vector<RingElement> res;
  RingElement tmpRingEl(ring);
  vector<Rcf> sv;
  for (int i = 0; i < (int)m; i++) {
    sv.assign(ringCoefficients.begin + i * n, ringCoefficients.begin + (i + 1) * n);
    tmpRingEl.setValues(sv, false); // Set them in the frequency domain
    res.push_back(tmpRingEl);
  }

  return res;
}
