///
///	\file   trapdoor.h
///
///	\brief  Class containing trapdoor generation and preimage sampling
///         algorithms for Ring LWE
///

#include <vector>

/// \class	Trapdoor
/// \brief	Class that implements the trapdoor algorithms from
///         "Implementation and Evaluation of Improved Gaussian Sampling for Lattice Trapdoors"
///          by K. Gur, Y. Polyakov, K. Rohloff, G. Ryan, and E. Savas

class Trapdoor : public ObjectBase {
public:
  Trapdoor(RNG *rng, std::shared_ptr<Ring> ring);
  ~Trapdoor;

  StatusCode TrapGen(double sigma, std::vector<RingElement> &A, std::vector<RingElement> &T);
  StatusCode GaussianPreimageSampling(std::vector<RingElement> &A, std::vector<RingElement> &T,
                                      RingElement &u, double sigma, std::vector<RingElement> &X);

  std::vector<RingElement> SampleG(double s, RingElement u);
  std::vector<RingElement> SamplePz(double s, double alpha, std::vector<RingElement> &T);

  Dim sizeOfA;
  Dim sizeOfT;
  double sampleStdev(double sigma);

  RNG *getRNG {
    return this->m_RNG_;
  }

protected:
  RNG *m_RNG_;
  std::shared_ptr<Ring> ring;
};
