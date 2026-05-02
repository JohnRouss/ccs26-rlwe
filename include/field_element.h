///
///	\file   field_element.h
///
///	\brief  Class to specify field elements. Each field element
///         is only specified by the dimension n.
///

#include <ostream>
#include <vector>

/// \class	FieldElement
/// \brief	Generic container for field LWE elements

class FieldElement : public ObjectBase, public std::vector<Fcf> {
public:
  FieldElement(Dim n);
  FieldElement(const FieldElement &elem) = default;
  FieldElement(FieldElement &&elem) noexcept = default;
  ~FieldElement = default;

  FieldElement &operator=(const FieldElement &rhs) = default;
  FieldElement &operator=(FieldElement &&rhs) noexcept = default;
  FieldElement &operator=(const Fcf &rhs);

  FieldElement &operator+=(const FieldElement &rhs);
  FieldElement &operator-=(const FieldElement &rhs);
  FieldElement &operator*=(const FieldElement &rhs); // Inefficient: Uses negacyclic convolution
  FieldElement &operator*=(const Fcf &scaleF); // Scaling all coefficients by scaleF factor (Faster
                                               // than multiplying by a constant field element)

  friend inline FieldElement operator+(FieldElement lhs, const FieldElement &rhs) {
    lhs += rhs;
    return lhs;
  };

  friend inline FieldElement operator-(FieldElement lhs, const FieldElement &rhs) {
    lhs -= rhs;
    return lhs;
  };

  friend inline FieldElement operator*(FieldElement lhs, const FieldElement &rhs) {
    lhs *= rhs;
    return lhs;
  };

  friend inline FieldElement operator*(FieldElement lhs, const Fcf &scaleF) {
    lhs *= scaleF;
    return lhs;
  };

  friend std::ostream &operator<<(std::ostream &os, FieldElement &elem);

  // Note: Use the following with care as it depends on LWE_EPS. It is not real equality
  friend bool operator==(const FieldElement &lhs, const FieldElement &rhs);
};

/// \class	FieldSpectrum
/// \brief	Generic container for field LWE spectra

class FieldSpectrum : public ObjectBase, public std::vector<FScf> {
public:
  FieldSpectrum(Dim n);
  FieldSpectrum(const FieldSpectrum &elem) = default;
  FieldSpectrum(FieldSpectrum &&elem) noexcept = default;
  ~FieldSpectrum = default;

  FieldSpectrum &operator=(const FieldSpectrum &rhs) = default;
  FieldSpectrum &operator=(FieldSpectrum &&rhs) noexcept = default;

  FieldSpectrum &operator+=(const FieldSpectrum &rhs);
  FieldSpectrum &operator-=(const FieldSpectrum &rhs);
  FieldSpectrum &operator*=(const FieldSpectrum &rhs);

  friend inline FieldSpectrum operator+(FieldSpectrum lhs, const FieldSpectrum &rhs) {
    lhs += rhs;
    return lhs;
  };

  friend inline FieldSpectrum operator-(FieldSpectrum lhs, const FieldSpectrum &rhs) {
    lhs -= rhs;
    return lhs;
  };

  friend inline FieldSpectrum operator*(FieldSpectrum lhs, const FieldSpectrum &rhs) {
    lhs *= rhs;
    return lhs;
  };

  FieldSpectrum operator~; // Inverse
  FieldSpectrum operator!; // Transposition

  friend std::ostream &operator<<(std::ostream &os, FieldSpectrum &elem);
};

// Transformation functions from field elements to their spectra
// Note: They take as an optional argument precomputed roots of unity
FieldSpectrum Fft(FieldElement &a, std::shared_ptr<FieldSpectrum> preCompZ = nullptr);
FieldElement Ifft(FieldSpectrum &a, std::shared_ptr<FieldSpectrum> preCompZ = nullptr);

// Other Helper functions
FieldElement ConvertIntegersToField(std::vector<Rcf> q);
void ExtractEvenAndOdd(const FieldElement &f, FieldElement &f0, FieldElement &f1);

