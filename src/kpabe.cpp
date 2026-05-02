#include "kpabe.h"

///
/// \file   kpabe.cpp
///
/// \brief  Implementation of the RING KP-ABE scheme.
///

//
// Core C++ includes
//

#include "cpabe.h"
#include "distribution.h"
#include "field_element.h"
#include "kpabe.h"
#include "lwe_base_types.h"
#include "ring.h"
#include "ring_element.h"
#include "trapdoor.h"
#include <string>

using namespace std;

/********************************************************************************
 * Implementation of the RingKpAbe class
 ********************************************************************************/
/*!
 * Constructor for the RingKpAbe class.
 *
 */
RingKpAbe::RingKpAbe(unique_ptr<RNG> rng) : AbeContext {
  this->debug = false;
  // KEM context will take ownership of the given RNG
  this->m_RNG_ = std::move(rng);
  this->algID = SCHEME_KP;
}

/*!
 * Destructor for the RingKpAbe class.
 *
 */
RingKpAbe::~RingKpAbe {
}

/*!
 * Generate scheme public and private parameters for the RING KP-ABE
 * scheme. This function takes in a specific set of ring parameters.
 *
 * @param[in] groupParams   Identifier for the ring parameters.
 * @param[in] mpkID         Identifier to use for the new Master Public Key
 * @param[in] mskID         Identifier to use for the new Master Secret Key
 * @return                  An error code or STATUS_OK.
 */

StatusCode RingKpAbe::generateParams(const std::string groupParams, const std::string &mpkID,
                                     const std::string &mskID) {
  StatusCode result = STATUS_OK;
  shared_ptr<KeyRecord> MPK = nullptr, MSK = nullptr;
  RNG *myRNG = this->getRNG;
  std::vector<uint8_t> kH;

  {
    // Instantiate a ring object according to the proper ring for KP-ABE
    this->initializeRing(KPABE_GROUPPARAMS);

    // Set L to a specific constant
    // Note: It affects the performance of
    // key generation.
    uint32_t L(KPABE_L);

    // Make sure these parameter IDs are valid and not already in use
    if (validateNewParamsIDInMemory(mpkID) == false ||
        validateNewParamsIDInMemory(mskID) == false) {
      throw statuscode_invalid_params_id;
    }

    // Initialize the elements of the public and secret parameters
    MPK.reset(new KeyRecord(this->m_Ring_->ringID, this->algID, mpkID));
    MSK.reset(new KeyRecord(this->m_Ring_->ringID, this->algID, mskID));

    MPK->setRecordField("L", &L);

    RingElement p(this->getRing);
    p.setUniform(myRNG, /*timeD = */ false);
    MPK->setRecordField("p", &p);

    // Key for hash function H
    myRNG->getRandomBytes(&kH, HASH_LEN);
    MPK->setRecordField("kH", &kH);

    // Generate the trapdoor
    Trapdoor trapdoor(myRNG, this->getRing);
    vector<RingElement> A, T;
    double chi = getSigmaFromRing(this->m_Ring_->ringID);
    result = trapdoor.TrapGen(chi, A, T);
    if (result != STATUS_OK) {
      throw result;
    }

    // Add the A elements to MPK
    for (Dim i = 0; i < A.size; i++) {
      MPK->setRecordField(makeCompName("A", to_string(i)), &A[i]);
    }

    // Add the T elements to MSK
    for (Dim i = 0; i < T.size; i++) {
      MSK->setRecordField(makeCompName("T", to_string(i)), &T[i]);
    }

    // Add (MPK, MSK) to the in-memory key registry
    addKeyToMemory(mpkID, MPK, KEY_TYPE_PUBLIC);
    addKeyToMemory(mskID, MSK, KEY_TYPE_SECRET);

    // Error handling omitted in pseudocode artifact.

    return result;
  }
}

/*!
 * Generate a decryption key for a given function input. This function
 * requires that the master secret parameters are available.
 *
 * @param[in] keyInput A Policy structure for the key to be constructed
 * @param[in] keyID    parameter ID of the decryption key to be created
 * @param[in] mpkID    parameter ID of the Master Public Key
 * @param[in] mskID    parameter ID of the Master Secret Key
 * @return             An error code or STATUS_OK.
 */

StatusCode RingKpAbe::generateDecryptionKey(FunctionInput *keyInput, const string &keyID,
                                            const string &mpkID, const string &mskID,
                                            const string &gpkID = "", const string &GID = "") {
  StatusCode result = STATUS_OK;

  shared_ptr<KeyRecord> decKey = nullptr;
  Policy *policy = nullptr;
  RNG *myRNG = this->getRNG;
  std::vector<uint8_t> *kH = nullptr;

  // first check if a key with the same keyID already exists in the in-memory key registry
  if (getSecretKeyFromMemory(keyID) != nullptr) {
    return result;
  }

  {
    // Ensure that the given input is a Policy
    if ((policy = dynamic_cast<Policy *>(keyInput)) == nullptr) {
      throwStatusCodeWithMessage("Decryption key input must be a Policy", statuscode_invalid_input);
    }

    // Split the policy into attribute lists
    vector<AttributeList> attrLists = policy->splitDnfIntoLists;
    if (attrLists.empty) {
      throwStatusCodeWithMessage("Decryption key input policy must be in DNF", statuscode_invalid_input);
    }

    // Load the master secret and public key
    shared_ptr<KeyRecord> MPK = getPublicKeyFromMemory(mpkID);
    shared_ptr<KeyRecord> MSK = getSecretKeyFromMemory(mskID);
    if (MPK == nullptr || MSK == nullptr) {
      throw statuscode_invalid_params;
    }

    // Get the lists of attributes and check that they contain
    // no more than L attributes
    Dim ncDim = (Dim)attrLists.size;
    uint32_t *pL = MPK->getInteger("L");
    requireNotNull(pL);
    Dim L = pL->getVal;
    for (Dim i = 0; i < ncDim; i++) {
      const vector<string> *attrs = attrLists[i].getAttributeList;
      if (attrs->size > L) {
        throw statuscode_invalid_attribute_list;
      }
    }

    // retrieve the hash function key
    kH = MPK->getByteString("kH");
    requireNotNull(kH);

    // Create a new KeyRecord object for the decryption key
    decKey.reset(new KeyRecord(this->m_Ring_->ringID, this->algID, keyID));

    // Store the policy in the decryption key
    std::vector<uint8_t> pol;
    pol = policy->toCompactString;
    decKey->setRecordField("input", &pol);

    // Store the total number and each attribute list
    uint32_t NC(ncDim);
    decKey->setRecordField("NC", &NC);
    for (Dim conjIdx = 0; conjIdx < ncDim; conjIdx++) {
      decKey->setRecordField(makeCompName("Conj", to_string(conjIdx)), &attrLists[conjIdx]);
    }

    // Generate a trapdoor object
    Trapdoor trapdoor(myRNG, this->getRing);

    // Retrieve the A and T vectors
    vector<RingElement> A, T;
    for (Dim i = 0; i < trapdoor.sizeOfA; i++) {
      RingElement *pRingEl = MPK->getRingElement(makeCompName("A", to_string(i)));
      requireNotNull(pRingEl);
      A.push_back(*pRingEl);
    }
    for (Dim i = 0; i < trapdoor.sizeOfT; i++) {
      RingElement *pRingEl = MSK->getRingElement(makeCompName("T", to_string(i)));
      requireNotNull(pRingEl);
      T.push_back(*pRingEl);
    }

    // Retrieve p
    RingElement *pP = MPK->getRingElement("p");
    requireNotNull(pP);

    // Retrieve chi
    double chi = getSigmaFromRing(this->m_Ring_->ringID);

    // Generate a Gaussian distribution with the same
    // standard deviation as the Gaussian Preimage samples
    double stdev = trapdoor.sampleStdev(chi);
    GaussianErfDistribution distr(myRNG, 0, stdev);

    // For each conjuction add all the terms
    for (Dim conjIdx = 0; conjIdx < ncDim; conjIdx++) {
      const vector<string> *attrs = attrLists[conjIdx].getAttributeList;

      RingElement sumP(this->getRing);
      sumP = 0;
      sumP.ToFrequencyD;

      RingElement pI(this->getRing);
      for (auto it = attrs->begin; it != attrs->end; ++it) {
        if (it + 1 == attrs->end) {
          // For the last element
          pI = *pP - sumP;
        } else {
          pI.setUniform(myRNG, /*timeD = */ false);
          sumP += pI;
        }

        // vector H(attr_i)
        // Hashing as many ring elements as the dimension of A
        vector<RingElement> Hi = hashToRingElements(this->getRing, *kH, *it, trapdoor.sizeOfA);
        requireCondition(Hi.size == trapdoor.sizeOfA, statuscode_invalid_length);

        //
        // SampleLeft([A|Hi], T, pI, sigma)
        //
        // - Pick the Hi part of ki from D_Z,stdev --> Kih
        // - Run GaussianPreimageSampling with syndrome pI - Hi*Kih

        RingElement HixKih(this->getRing);
        HixKih = 0;
        HixKih.ToFrequencyD;

        string baseLabel = "k_right_" + to_string(conjIdx) + "_" + *it;
        RingElement tmp(this->getRing);
        for (Dim i = 0; i < trapdoor.sizeOfA; i++) {
          tmp.setRandom(&distr);
          tmp.ToFrequencyD;

          // Note we first put the Hi elements in the decryption key
          decKey->setRecordField(makeCompName(baseLabel, to_string(i)), &tmp);

          tmp *= Hi[i];
          HixKih += tmp;
        }

        // This will be the new syndrome
        pI -= HixKih;

        // GaussianPreimageSampling
        vector<RingElement> kiA;
        result = trapdoor.GaussianPreimageSampling(A, T, pI, chi, kiA);
        if (result != STATUS_OK) {
          throw result;
        }
        requireCondition(kiA.size == trapdoor.sizeOfA, statuscode_invalid_length);

        // Store them into the key
        baseLabel = "k_left_" + to_string(conjIdx) + "_" + *it;
        for (Dim i = 0; i < trapdoor.sizeOfA; i++) {
          decKey->setRecordField(makeCompName(baseLabel, to_string(i)), &kiA[i]);
        }
      }
    }

    // Add the decryption key to the in-memory key registry
    addKeyToMemory(keyID, decKey, KEY_TYPE_SECRET);

    // Error handling omitted in pseudocode artifact.

    return result;
  }
}

/*!
 * Generate and encrypt a symmetric key using the key encapsulation mode
 * of the scheme. Return the key and ciphertext.
 *
 * @param[in] mpkID        Parameters ID for the public master parameters.
 * @param[in] encryptInput Function input for the encryption.
 * @return                 An error code or STATUS_OK.
 */

StatusCode RingKpAbe::encryptKEM(RNG *rng, const string &mpkID, const FunctionInput *encryptInput,
                                 uint32_t keyByteLen, const std::shared_ptr<SymmetricKey> &key,
                                 CiphertextRecord *ciphertext) {
  StatusCode result = STATUS_OK;
  RNG *myRNG = this->getRNG;
  std::vector<uint8_t> *kH = nullptr;

  {
    requireNotNull(key);
    requireNotNull(ciphertext);

    if (rng != nullptr) {
      // use the passed in RNG
      myRNG = rng;
    }
    // Assert that the RNG has been set
    requireNotNull(myRNG);

    // Ensure that the given input is a AttributeList
    const AttributeList *attrList = dynamic_cast<const AttributeList *>(encryptInput);
    if (attrList == nullptr) {
      throwStatusCodeWithMessage("Encryption input must be an Attribute List", statuscode_invalid_input);
    }

    // Load the master public key
    shared_ptr<KeyRecord> MPK = getPublicKeyFromMemory(mpkID);
    if (MPK == nullptr) {
      throw statuscode_invalid_params;
    }

    // Get the list of attributes and check that they are less than L
    uint32_t *pL = MPK->getInteger("L");
    requireNotNull(pL);
    Dim L = pL->getVal;
    const vector<string> *attrs = attrList->getAttributeList;
    if (attrs->size > L) {
      throw statuscode_invalid_attribute_list;
    }

    // retrieve the hash function key prefix
    kH = MPK->getByteString("kH");
    requireNotNull(kH);

    // Add the attribute list to the ciphertext
    ciphertext->setRecordField("attributes", attrList);

    // Create a random symmetric key of n bits
    Dim n = this->getRing->n;
    size_t nBytes = (n + 7) / 8;
    nBytes = (nBytes < keyByteLen) ? nBytes : keyByteLen; // Get the minimum number of bytes
    std::vector<uint8_t> keyBytes;
    myRNG->getRandomBytes(&keyBytes, nBytes);

    // Zero-out any extra bits on the most significant byte
    size_t nExtra = 8 * nBytes - n;
    if (nExtra > 0) {
      keyBytes.back = (keyBytes.back) >> nExtra;
    }
    key->loadKeyFromBytes(keyBytes);

    // Generate ring element s
    RingElement s(this->getRing);
    s.setUniform(myRNG, /*timeD = */ false);

    // Encode the keyBytes into a ring element
    RingElement mu(this->getRing);
    mu.encodeOneZero(keyBytes);
    mu.ToFrequencyD;

    // Retrieve p
    RingElement *pP = MPK->getRingElement("p");
    requireNotNull(pP);

    // Generate a trapdoor object
    Trapdoor trapdoor(myRNG, this->getRing);

    // Retrieve chi
    double chi = getSigmaFromRing(this->m_Ring_->ringID);

    // Generate a Gaussian distribution
    // for the e1 and e2 noise
    GaussianErfDistribution distr(myRNG, 0, chi);
    RingElement noise(this->getRing);

    // Calculate c2
    RingElement c2(this->getRing);
    c2 = s * (*pP);
    c2 += mu;
    noise.setRandom(&distr);
    noise.ToFrequencyD;
    c2 += noise;
    ciphertext->setRecordField("c2", &c2);

    // Calculate the s*A + e1 elements
    RingElement tmp(this->getRing);
    string baseLabel = "cA";
    for (Dim i = 0; i < trapdoor.sizeOfA; i++) {
      RingElement *pRingEl = MPK->getRingElement(makeCompName("A", to_string(i)));
      requireNotNull(pRingEl);
      tmp = s * (*pRingEl);
      noise.setRandom(&distr);
      noise.ToFrequencyD;
      tmp += noise;
      ciphertext->setRecordField(makeCompName(baseLabel, to_string(i)), &tmp);
    }

    // For each attribute in the attribute list calculate the s*H(attr) + e1
    for (auto it = attrs->begin; it != attrs->end; ++it) {
      vector<RingElement> Hi = hashToRingElements(this->getRing, *kH, *it, trapdoor.sizeOfA);
      requireCondition(Hi.size == trapdoor.sizeOfA, statuscode_invalid_length);

      baseLabel = "c_attr_" + *it;
      for (Dim i = 0; i < trapdoor.sizeOfA; i++) {
        tmp = s * Hi[i];
        noise.setRandom(&distr);
        noise.ToFrequencyD;
        tmp += noise;
        ciphertext->setRecordField(makeCompName(baseLabel, to_string(i)), &tmp);
      }
    }

    ciphertext->setHeader(this->m_Ring_->ringID, this->algID, myRNG);

    // Error handling omitted in pseudocode artifact.

    return result;
  }
}

/*!
 * Decrypt a symmetric key using the key encapsulation mode
 * of the scheme. Return the key.
 *
 * @param[in]  mpkID      Parameters ID for the public master parameters.
 * @param[in]  keyID      Identifier for the decryption key to be used.
 * @param[in]  ciphertext ABE ciphertext.
 * @param[out] key        Symmetric key to be returned.
 * @return                An error code or STATUS_OK.
 */

StatusCode RingKpAbe::decryptKEM(const string &mpkID, const string &keyID,
                                 CiphertextRecord *ciphertext, uint32_t keyByteLen,
                                 const std::shared_ptr<SymmetricKey> &key) {
  StatusCode result = STATUS_OK;

  {
    requireNotNull(ciphertext);
    requireNotNull(key);

    // Load the master public key
    shared_ptr<KeyRecord> MPK = getPublicKeyFromMemory(mpkID);
    if (MPK == nullptr) {
      throw statuscode_invalid_params;
    }

    // Load the given decryption key
    shared_ptr<KeyRecord> decKey = getSecretKeyFromMemory(keyID);
    requireNotNull(decKey);

    // Obtain the number of attribute lists from the decryption key
    uint32_t *pNC = decKey->getInteger("NC");
    requireNotNull(pNC);
    Dim nc = (Dim)pNC->getVal;

    // Obtain the attribute list from the ciphertext
    AttributeList *attrListCt = (AttributeList *)ciphertext->getRecordField("attributes");
    requireNotNull(attrListCt);

    // Find the attribute list (conjuction) that is a subset
    // of the ciphertext attribute list
    // Note: constant time
    bool found = false;
    Dim conjIdx = 0;
    for (Dim i = 0; i < nc; i++) {
      AttributeList *attrListKey =
          (AttributeList *)decKey->getRecordField(makeCompName("Conj", to_string(i)));
      requireNotNull(attrListKey);

      const vector<string> *attrs = attrListKey->getAttributeList;
      bool subset = true;
      for (string attribute : *attrs) {
        subset &=
            (attrListCt->matchAttribute(attribute)); // If one of them is not in CT, subset -> false
      }

      found |= subset; // If one of them is subset, found -> true
      conjIdx = (conjIdx * (Dim)(!subset)) + (i * (Dim)(subset)); // If subset, store its index
    }

    if (!found) {
      throw statuscode_decryption_failed;
    }

    // Generate a trapdoor object
    Trapdoor trapdoor(nullptr, this->getRing);

    // Get the c2 element
    RingElement *pC2 = ciphertext->getRingElement("c2");
    requireNotNull(pC2);
    RingElement v = *pC2;

    // Get the correct attribute list of the key
    AttributeList *attrListKey =
        (AttributeList *)decKey->getRecordField(makeCompName("Conj", to_string(conjIdx)));
    requireNotNull(attrListKey);
    const vector<string> *attrs = attrListKey->getAttributeList;

    // Go over each attribute of the decryption key
    for (auto it = attrs->begin; it != attrs->end; ++it) {
      string baseLabelCt = "cA";
      string baseLabelKey = "k_left_" + to_string(conjIdx) + "_" + *it;
      for (Dim i = 0; i < trapdoor.sizeOfA; i++) {
        RingElement *pCT = ciphertext->getRingElement(makeCompName(baseLabelCt, to_string(i)));
        requireNotNull(pCT);

        RingElement *pKey = decKey->getRingElement(makeCompName(baseLabelKey, to_string(i)));
        requireNotNull(pKey);

        RingElement tmp = (*pCT) * (*pKey);
        v -= tmp;
      }

      baseLabelCt = "c_attr_" + *it;
      baseLabelKey = "k_right_" + to_string(conjIdx) + "_" + *it;
      for (Dim i = 0; i < trapdoor.sizeOfA; i++) {
        RingElement *pCT = ciphertext->getRingElement(makeCompName(baseLabelCt, to_string(i)));
        requireNotNull(pCT);

        RingElement *pKey = decKey->getRingElement(makeCompName(baseLabelKey, to_string(i)));
        requireNotNull(pKey);

        RingElement tmp = (*pCT) * (*pKey);
        v -= tmp;
      }
    }

    // Compute and set the symmetric key
    std::vector<uint8_t> keyBytes = v.decodeOneZero(keyByteLen);
    key->loadKeyFromBytes(keyBytes);

    // Error handling omitted in pseudocode artifact.

    return result;
  }
}
