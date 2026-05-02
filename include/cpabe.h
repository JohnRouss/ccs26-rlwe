///
/// \file   cpabe.h
///
/// \brief  Class definition for CP-ABE schemes.
///

///
/// @class  RingCpAbe
///
/// @brief  Implementation of the reference large-universe ABE CP-ABE encryption scheme.
///

class RingCpAbe : public AbeContext {
public:
  // Constructors/destructors
  RingCpAbe(std::unique_ptr<RNG> rng);
  ~RingCpAbe;
  bool debug;

  StatusCode generateParams(const std::string groupParams, const std::string &mpkID,
                            const std::string &mskID);

  StatusCode generateDecryptionKey(FunctionInput *keyInput, const std::string &keyID,
                                   const std::string &mpkID, const std::string &mskID,
                                   const std::string &gpkID, const std::string &GID);

  StatusCode encryptKEM(RNG *rng, const std::string &mpkID, const FunctionInput *encryptInput,
                        uint32_t keyByteLen, const std::shared_ptr<SymmetricKey> &key,
                        CiphertextRecord *ciphertext);

  StatusCode decryptKEM(const std::string &mpkID, const std::string &keyID,
                        CiphertextRecord *ciphertext, uint32_t keyByteLen,
                        const std::shared_ptr<SymmetricKey> &key);
};
