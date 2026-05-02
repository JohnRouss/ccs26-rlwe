#include "ring.h"

///
///	\file   ring.cpp
///
///	\brief  Class implementation for rings.
///

#include "cpabe.h"
#include "distribution.h"
#include "field_element.h"
#include "kpabe.h"
#include "lwe_base_types.h"
#include "ring.h"
#include "ring_element.h"
#include "trapdoor.h"

using namespace std;

// Utility functions

// Get the ring ID from the ringParams string
// Note: If the default BP params are specified
//       pick the R_60 ring.
RingPresetId getRingID(const string &ringParams) {
  RingPresetId curveID = RING_PRESET_NONE;

  if (ringParams == KPABE_GROUPPARAMS) {
    curveID = RING_PRESET_PRODUCTION_57;
  } else if (ringParams == CPABE_GROUPPARAMS) {
    curveID = RING_PRESET_PRODUCTION_60;
  } else {
    // Unrecognized parameter type
    throw statuscode_invalid_params;
  }

  return curveID;
}

// Get the ring parameters structure from the
// ring ID
RingParameters getRingParameters(RingPresetId ringID) {
  RingParameters res; // Initialize with default 0 values
  for (auto params : PresetRings::Parameters) {
    if (params.ringID == ringID) {
      res = params;
      break;
    }
  }
  return res;
}

// Global function for creating a new ring object
// according to the preset rings
Ring *createNewRing(const string ringParams) {
  RingPresetId ringID = getRingID(ringParams);

  RingParameters ringParamsStruct = getRingParameters(ringID);
  if (!ringParamsStruct.production) {
    return nullptr;
  }

  Ring *res = nullptr;
  res = new Ring(ringParamsStruct);
  res->ringID = ringID;

  return res;
}

// Global function to get the appropriate sigma
double getSigmaFromRing(RingPresetId ringID) {
  double res = -1.0; // Invalid result
  for (auto params : PresetRings::Parameters) {
    if (params.ringID == ringID) {
      res = params.sigma;
      break;
    }
  }
  if (res < 0.0) {
    throw statuscode_invalid_params;
  } else {
    return res;
  }
}

// Helper function to invert ring coefficients
// mod q via the Extended Euclidean algorithm
Rcf RcfInvertModQ(Rcf x, Rcf q) {
  if (x >= q) {
    throw statuscode_invalid_params;
    return 0;
  }

  Rcf r0 = q;
  Rcf r1 = x;
  Rcf t0 = 0;
  Rcf t1 = 1;

  while (r1 > 1) {
    Rcf rn = r0 % r1; // Remainder
    Rcf qn = r0 / r1; // Quotient

    // tn = t0 - qn * t1 mod q
    Rcf prod = RMul(qn, t1, q);
    Rcf tn = RSub(t0, prod, q);

    // Move one step
    r0 = r1;
    t0 = t1;
    r1 = rn;
    t1 = tn;
  }

  if (r1 == 0) {
    throw statuscode_divide_by_zero;
    return 0;
  }

  return t1;
}

// Number of bits of a number
Dim log2D(Rcf value) {
  Dim res = 0;
  while (value > 0) {
    value >>= 1;
    res++;
  }
  return res;
}

// Check if a dimension is a power of two
void assertPowerOf2(Dim n) {
  if ((n == 0) || ((n & (n - 1)) != 0)) {
    throw statuscode_invalid_params;
  }
}

/********************************************************************************
 * Implementation of the Ring class
 ********************************************************************************/
/*!
 * Constructors for the Ring class.
 *
 */

Ring::Ring(Rcf q, Dim n) : q(q), n(n) {
  this->ringID = RING_PRESET_NONE;

  assertPowerOf2(n);

  this->logn = log2D(n - 1);
  this->nInv = RcfInvertModQ((Rcf)n, q);
  this->logq = log2D(q);
};

Ring::Ring(Rcf q, Dim n, Rcf psi) : Ring(q, n) {
  this->setPsi(psi);
};

Ring::Ring(RingParameters p) {
  this->ringID = p.ringID;

  assertPowerOf2(p.n);

  this->n = p.n;
  this->logn = log2D(p.n - 1);

  Rcf Q = 0;
  Rcf Psi = 0;

  // Fix q and psi from the two 64-bit literals
#ifndef __RCF64__
  if (p.qHi > 0) {
    Q = (Rcf)p.qHi;
    Q <<= 64;
  }
#endif
  Q |= (Rcf)p.qLo;

#ifndef __RCF64__
  if (p.psiHi > 0) {
    Psi = (Rcf)p.psiHi;
    Psi <<= 64;
  }
#endif
  Psi |= (Rcf)p.psiLo;

  this->q = Q;
  this->nInv = RcfInvertModQ((Rcf)p.n, Q);
  this->logq = log2D(Q);

  this->setPsi(Psi);
};

Ring::Ring(RingPresetId ringID) : Ring(getRingParameters(ringID)) {};

/*!
 * Destructor for the Ring class.
 *
 */

Ring::~Ring {
}

void Ring::setPsi(Rcf psi) {
  this->psi = psi;
  this->psiInv = RcfInvertModQ(psi, q);

  // Store all the powers of Psi (up to n-1)
  // in normal order
  vector<Rcf> tmpPsiVec;
  Rcf tmp = 1;
  tmpPsiVec.push_back(tmp);
  for (Dim i = 0; i < this->n - 1; i++) {
    tmp = RMul(tmp, psi, q);
    tmpPsiVec.push_back(tmp);
  }

  // And the inverse powers of Psi
  vector<Rcf> tmpPsiInvVec;
  Rcf psiInv = this->psiInv;
  tmp = 1;
  tmpPsiInvVec.push_back(tmp);
  for (Dim i = 0; i < this->n - 1; i++) {
    tmp = RMul(tmp, psiInv, q);
    tmpPsiInvVec.push_back(tmp);
  }

  // Copy them in bit-reversed order
  this->psiVecR.clear;
  this->psiInvVecR.clear;
  for (Dim i = 0; i < this->n; i++) {
    Rcf idx = BRev(i, this->logn);
    this->psiVecR.push_back(tmpPsiVec[idx]);
    this->psiInvVecR.push_back(tmpPsiInvVec[idx]);
  }

  tmpPsiVec.clear;
  tmpPsiInvVec.clear;
}

ostream &operator<<(ostream &os, Ring &ring) {
  os << "(q, logq, n, logn, n^{-1}";
  if (ring.psi != 0) {
    os << ", psi, psi^{-1}";
  }
  os << ") = (";

  os << ring.q << ", " << ring.logq << ", " << ring.n << ", " << ring.logn << ", " << ring.nInv;
  if (ring.psi != 0) {
    os << ", " << ring.psi << ", " << ring.psiInv;
  }
  os << ")";

  return os;
}
