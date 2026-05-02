#include "field_element.h"

///
///	\file   field_element.cpp
///
///	\brief  Class implementation for field elements.
///

#include "cpabe.h"
#include "distribution.h"
#include "field_element.h"
#include "kpabe.h"
#include "lwe_base_types.h"
#include "ring.h"
#include "ring_element.h"
#include "trapdoor.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <ostream>
#include <vector>

using namespace std;

/********************************************************************************
 * Implementation of the FieldElement class
 ********************************************************************************/
/*!
 * Constructor for the FieldElement class.
 *
 */
FieldElement::FieldElement(Dim n) : std::vector<Fcf>(n, 0) {
}

FieldElement &FieldElement::operator=(const Fcf &rhs) {
  Dim n = (Dim)this->size;

  (*this)[0] = rhs;
  for (Dim i = 1; i < n; i++) {
    (*this)[i] = 0;
  }

  return *this;
}

FieldElement &FieldElement::operator+=(const FieldElement &rhs) {
  if (this->size != rhs.size) {
    throw statuscode_invalid_field_element_size;
  }

  const Dim n = (Dim)this->size;
  for (Dim i = 0; i < n; i++) {
    (*this)[i] += rhs[i];
  }

  return *this;
}

FieldElement &FieldElement::operator-=(const FieldElement &rhs) {
  if (this->size != rhs.size) {
    throw statuscode_invalid_field_element_size;
  }

  const Dim n = (Dim)this->size;
  for (Dim i = 0; i < n; i++) {
    (*this)[i] -= rhs[i];
  }

  return *this;
}

FieldElement &FieldElement::operator*=(const FieldElement &rhs) {
  if (this->size != rhs.size) {
    throw statuscode_invalid_field_element_size;
  }

  const Dim n = (Dim)this->size;
  FieldElement res(n);

  // First make a reverse copy of the rhs vector
  vector<Fcf> revRhs(rhs.rbegin, rhs.rend);

  // Then go through both vectors multiplying elements and adding (or subtracting)
  // monomials according to the negative wrapped convolution.
  for (Dim k = 0; k < n; k++) {
    Fcf acc = 0;

    // First put the last element of reverse rhs in front
    std::rotate(revRhs.begin, revRhs.end - 1, revRhs.end);

    for (Dim i = 0; i < n; i++) {
      Fcf monomial = (*this)[i] * revRhs[i];

      if (i > k) {
        // Subtract
        acc -= monomial;
      } else {
        // Add
        acc += monomial;
      }
    }

    res[k] = acc;
  }

  // Copy the vector
  this->assign(res.begin, res.end);

  return *this;
}

FieldElement &FieldElement::operator*=(const Fcf &scaleF) {
  const Dim n = (Dim)this->size;

  for (Dim i = 0; i < n; i++) {
    (*this)[i] = scaleF * (*this)[i];
  }

  return *this;
}

ostream &operator<<(ostream &os, FieldElement &elem) {
  os << "[ ";
  for (auto cf : elem) {
    os << cf << " ";
  }
  os << "]";
  return os;
}

bool operator==(const FieldElement &lhs, const FieldElement &rhs) {
  if (lhs.size != rhs.size) {
    return false;
  }

  bool rc = true;
  for (size_t i = 0; i < lhs.size; i++) {
    Fcf diff = std::fabs(lhs[i] - rhs[i]);
    rc &= (diff < LWE_EPS);
  }

  return rc;
}

/********************************************************************************
 * Implementation of the FieldSpectrum class
 ********************************************************************************/
/*!
 * Constructor for the FieldSpectrum class.
 *
 */
FieldSpectrum::FieldSpectrum(Dim n) : std::vector<FScf>(n, 0) {
}

FieldSpectrum &FieldSpectrum::operator+=(const FieldSpectrum &rhs) {
  if (this->size != rhs.size) {
    throw statuscode_invalid_field_element_size;
  }

  const Dim n = (Dim)this->size;
  for (Dim i = 0; i < n; i++) {
    (*this)[i] += rhs[i];
  }

  return *this;
}

FieldSpectrum &FieldSpectrum::operator-=(const FieldSpectrum &rhs) {
  if (this->size != rhs.size) {
    throw statuscode_invalid_field_element_size;
  }

  const Dim n = (Dim)this->size;
  for (Dim i = 0; i < n; i++) {
    (*this)[i] -= rhs[i];
  }

  return *this;
}

FieldSpectrum &FieldSpectrum::operator*=(const FieldSpectrum &rhs) {
  if (this->size != rhs.size) {
    throw statuscode_invalid_field_element_size;
  }

  const Dim n = (Dim)this->size;
  for (Dim i = 0; i < n; i++) {
    (*this)[i] *= rhs[i];
  }

  return *this;
}

FieldSpectrum FieldSpectrum::operator~{
  // Invert by calculating 1/(a+bi) for each
  // coordinate in the spectrum. This is
  // equal to (a-bi) / norm(a+bi).

  const Dim n = (Dim)this->size;
  FieldSpectrum res(n);
  for (Dim k = 0; k < n; k++) {
    Fcf denom = std::norm((*this)[k]);
    if (denom < LWE_EPS) {
      throw statuscode_divide_by_zero;
    }
    res[k] = std::conj((*this)[k]) / denom;
  }

  return res;
}

FieldSpectrum FieldSpectrum::operator!{
  // Transpose by permuting the indices of each
  // coordinate k -> n - 1 - k.
  // This is because tranposition in the canonical
  // embedding maps zeta^{2k + 1} to zeta^{-(2k+1)} =
  // zeta^{2(n - k - 1) + 1}

  const Dim n = (Dim)this->size;
  FieldSpectrum res(n);
  for (Dim k = 0; k < n; k++) {
    res[k] = (*this)[n - 1 - k];
  }

  return res;
}

ostream &operator<<(ostream &os, FieldSpectrum &elem) {
  os << "< ";
  for (auto cf : elem) {
    os << cf << " ";
  }
  os << ">";
  return os;
}

/********************************************************************************
 * FFT transformations
 ********************************************************************************/

// General function for FFT and inverse FFT
static inline void fft_inplace(FieldSpectrum &a, bool inverse) {
  const Dim n = (Dim)a.size;
  if ((n == 0) || ((n & (n - 1)) != 0)) {
    throw statuscode_invalid_length;
  }

  Dim logn = 0;
  Dim nTmp = n >> 1;
  while (nTmp > 0) {
    logn++;
    nTmp >>= 1;
  }

  // Bit-reverse the input vector
  for (Dim i = 1; i < n; i++) {
    Dim j = BRev(i, logn);
    if (i < j) {
      std::swap(a[i], a[j]);
    }
  }

  // FFT stages
  for (Dim len = 2; len <= n; len <<= 1) {
    const Fcf ang = (inverse ? -1.0 : 1.0) * (2.0 * M_PI / (Fcf)len);
    const complex<Fcf> wlen(cos(ang), sin(ang));
    for (Dim i = 0; i < n; i += len) {
      complex<Fcf> w(1.0, 0.0);
      const Dim half = len >> 1;
      for (Dim j = 0; j < half; ++j) {
        auto u = a[i + j];
        auto v = a[i + j + half] * w;
        a[i + j] = u + v;
        a[i + j + half] = u - v;
        w *= wlen;
      }
    }
  }

  // For the inverse FFT divide by n
  if (inverse) {
    const Fcf invn = 1.0 / (Fcf)n;
    for (auto &z : a)
      z *= invn;
  }
}

// Helper functions for phases zeta^m
static inline complex<Fcf> cis(Fcf theta) {
  return {std::cos(theta), std::sin(theta)};
}

// Precompute exp(i*2π/2n) once if you batch many calls.
static inline complex<Fcf> zeta_pow(Dim m, Dim n2) {
  const Fcf theta = (2.0 * M_PI / (Fcf)n2) * (Fcf)m;
  return cis(theta);
}

FieldSpectrum Fft(FieldElement &a, shared_ptr<FieldSpectrum> preCompZ) {
  Dim n = (Dim)a.size;
  Dim n2 = n << 1;

  FieldSpectrum buf(n);
  // Pre-twist by zeta^i, then FFT length d (omega = e^{2 Pi i / n})
  for (Dim i = 0; i < n; ++i) {
    complex<Fcf> phase;
    if (preCompZ == nullptr) {
      phase = zeta_pow(i, n2); // zeta^i = e^{i 2 Pi i / n2}
    } else {
      phase = preCompZ->at(i);
    }
    buf[i] = phase * a[i]; // complex
  }

  fft_inplace(buf, /*inverse=*/false);
  // buf[k] = sum_i a_i zeta^i (zeta^2)^{k i} = sum_i a_i zeta^{(2k+1)i}  (evaluation at
  // zeta^{2k+1})
  return buf; // size n
}

FieldElement Ifft(FieldSpectrum &a, shared_ptr<FieldSpectrum> preCompZ) {
  Dim n = (Dim)a.size;
  Dim n2 = n << 1;

  FieldSpectrum buf = a;
  fft_inplace(buf, /*inverse=*/true); // inverse FFT with pre-twisted time domain

  FieldElement res(n);
  for (Dim i = 0; i < n; ++i) {
    complex<Fcf> phase_inv;
    if (preCompZ == nullptr) {
      phase_inv = std::conj(zeta_pow(i, n2)); // zeta^{-i}
    } else {
      phase_inv = preCompZ->at(i);
    }

    auto val = buf[i] * phase_inv; // should be (real) a_i

    res[i] = val.real; // discard tiny imag noise
  }
  return res;
}

FieldElement ConvertIntegersToField(std::vector<Rcf> q) {
  Dim n = (Dim)q.size;
  FieldElement res(n);
  for (Dim k = 0; k < n; k++) {
    res[k] = (Fcf)q[k];
  }
  return res;
}

void ExtractEvenAndOdd(const FieldElement &f, FieldElement &f0, FieldElement &f1) {
  Dim n = (Dim)f.size;
  Dim n2_0 = (Dim)f0.size;
  Dim n2_1 = (Dim)f1.size;

  requireCondition(((n = 2 * n2_0) && (n2_0 == n2_1)), statuscode_invalid_length);

  for (Dim k = 0; k < n; k++) {
    if (k % 2 == 0) {
      f0[k / 2] = f[k];
    } else {
      f1[k / 2] = f[k];
    }
  }
}
