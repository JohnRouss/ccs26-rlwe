///
///	\file   ring_element.h
///
///	\brief  Class to specify ring elements. Each ring element
///         is linked to a Ring object.
///

#include <ostream>
#include <vector>

/// \class	RingElement
/// \brief	Generic container for ring LWE elements

class RingElement : public ObjectBase, public std::vector<Rcf> {
protected:
  std::shared_ptr<Ring> ring; // A pointer to the ring structure of the element
  bool timeD; // Whether it is on time or frequency domain (If false, then the coefficients are
              // stored in bit reversed order)

public:
  RingElement(std::shared_ptr<Ring> ring, bool timeDomain = true);
  RingElement(const RingElement &elem) = default;
  RingElement(RingElement &&elem) noexcept = default;
  ~RingElement = default;

  std::shared_ptr<Ring> getRing const;
  bool isTimeDomain const;

  void setValues(std::vector<Rcf> &input, bool timeDomain = true);
  void setConstantTerm(Rcf c, bool timeDomain = true);
  void setRandom(Distribution *distr, bool timeDomain = true);
  void setUniform(RNG *rng, bool timeDomain = true);

  // Encode a string of ones and zeros into a
  // ring element and the inverse
  void encodeOneZero(std::vector<uint8_t> &msg);
  std::vector<uint8_t> decodeOneZero(uint32_t keyByteLen);

  // Number theoretic transform and their inverse
  void NaiveNtt;
  void NaiveIntt;
  void FastNtt;
  void FastIntt;

  // Convert to specific domains
  void ToTimeD;
  void ToFrequencyD;

  // Operators
  RingElement &operator=(const RingElement &rhs) = default;
  RingElement &operator=(RingElement &&rhs) noexcept = default;
  RingElement &operator=(const Rcf &rhs);

  RingElement &operator+=(const RingElement &rhs);
  RingElement &operator-=(const RingElement &rhs);
  RingElement &operator*=(const RingElement &rhs);
  RingElement &operator*=(const Rcf &scaleR);

  friend inline RingElement operator+(RingElement lhs, const RingElement &rhs) {
    lhs += rhs;
    return lhs;
  };

  friend inline RingElement operator-(RingElement lhs, const RingElement &rhs) {
    lhs -= rhs;
    return lhs;
  };

  friend inline RingElement operator*(RingElement lhs, const RingElement &rhs) {
    lhs *= rhs;
    return lhs;
  };

  friend inline RingElement operator*(RingElement lhs, const Rcf &scaleR) {
    lhs *= scaleR;
    return lhs;
  };

  RingElement operator!; // Transposition

  friend std::ostream &operator<<(std::ostream &os, RingElement &elem);
  friend bool operator==(const RingElement &lhs, const RingElement &rhs);
  bool isEqual(ObjectBase *z) const;

  RingElement *clone const {
    return new RingElement(*this);
  }

  void serialize(std::vector<uint8_t> &result) const;

  // Return the infinity norm with center in 0
  Rcf InfNorm;
};

// Helper function to hash a string into a vector of m ring elements
std::vector<RingElement> hashToRingElements(std::shared_ptr<Ring> ring, std::vector<uint8_t> &key,
                                            const std::string &input, Dim m);
