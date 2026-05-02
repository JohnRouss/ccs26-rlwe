#include "cpabe.h"

///
/// \file   cpabe.cpp
///
/// \brief  Implementation of the RING CP-ABE scheme.
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
 * Implementation of the RingCpAbe class
 ********************************************************************************/
/*!
 * Constructor for the RingCpAbe class.
 *
 */
RingCpAbe::RingCpAbe(unique_ptr<RNG> rng) : AbeContext {
  this->debug = false;
  // KEM context will take ownership of the given RNG
  this->m_RNG_ = std::move(rng);
  this->algID = SCHEME_CP;
}

/*!
 * Destructor for the RingCpAbe class.
 *
 */
RingCpAbe::~RingCpAbe {
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

StatusCode RingCpAbe::generateParams(const std::string groupParams, const std::string &mpkID,
                                     const std::string &mskID) {
  StatusCode result = STATUS_OK;
  shared_ptr<KeyRecord> MPK = nullptr, MSK = nullptr;
  RNG *myRNG = this->getRNG;
  std::vector<uint8_t> kH;

  // Instantiate a ring object according to the proper ring for CP-ABE
  this->initializeRing(CPABE_GROUPPARAMS);

  // Set L to a specific constant
  // Note: It affects the performance of
  // key generation.
  uint32_t L(CPABE_L);

  // Make sure these parameter IDs are valid and not already in use
  if (validateNewParamsIDInMemory(mpkID) == false || validateNewParamsIDInMemory(mskID) == false) {
    throw statuscode_invalid_params_id;
  }

  // Initialize the elements of the public and secret parameters
  MPK.reset(new KeyRecord(this->m_Ring_->ringID, this->algID, mpkID));
  MSK.reset(new KeyRecord(this->m_Ring_->ringID, this->algID, mskID));

  MPK->setRecordField("L", &L);

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

StatusCode RingCpAbe::generateDecryptionKey(FunctionInput *keyInput, const string &keyID,
                                            const string &mpkID, const string &mskID,
                                            const string &gpkID = "", const string &GID = "") {
  StatusCode result = STATUS_OK;

  shared_ptr<KeyRecord> decKey = nullptr;
  RNG *myRNG = this->getRNG;
  std::vector<uint8_t> *kH = nullptr;

  // first check if a key with the same keyID already exists in the in-memory key registry
  if (getSecretKeyFromMemory(keyID) != nullptr) {
    return result;
  }

  // Ensure that the given input is a AttributeList
  const AttributeList *attrList = dynamic_cast<const AttributeList *>(keyInput);
  if (attrList == nullptr) {
    throwStatusCodeWithMessage("Decryption key input must be an Attribute List",
                          statuscode_invalid_input);
  }

  // Load the master secret and public key
  shared_ptr<KeyRecord> MPK = getPublicKeyFromMemory(mpkID);
  shared_ptr<KeyRecord> MSK = getSecretKeyFromMemory(mskID);
  if (MPK == nullptr || MSK == nullptr) {
    throw statuscode_invalid_params;
  }

  // Get the list of attributes and check that they are fewer than L
  uint32_t *pL = MPK->getInteger("L");
  requireNotNull(pL);
  Dim L = pL->getVal;
  const vector<string> *attrs = attrList->getAttributeList;
  Dim t = (Dim)(attrs->size);
  if (t > L) {
    throw statuscode_invalid_attribute_list;
  }

  // retrieve the hash function key
  kH = MPK->getByteString("kH");
  requireNotNull(kH);

  // Create a new KeyRecord object for the decryption key
  decKey.reset(new KeyRecord(this->m_Ring_->ringID, this->algID, keyID));

  // Store the attribute list in the decryption key
  decKey->setRecordField("input", attrList);

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

  // Retrieve chi
  double chi = getSigmaFromRing(this->m_Ring_->ringID);

  // Generate a Gaussian distribution with chi standard deviation
  GaussianErfDistribution distrR(myRNG, 0, chi);

  // Generate the all the r_ell elements
  // Total number: m + 3m + 5m + ... +(2L-1)m =
  // m * (2L) * L / 2 = m*L^2
  // And fill a bytestring with their bytes
  vector<Rcf> rell;
  Dim m = trapdoor.sizeOfA;
  Dim numOfRell = m * L * L;
  vector<uint8_t> rellBytes;
  for (Dim i = 0; i < numOfRell; i++) {
    Rcf tmpRell = distrR.Sample;
    rell.push_back(tmpRell);

    size_t numOfBytes = sizeof(Rcf);
    while (numOfBytes-- > 0) {
      rellBytes.push_back((uint8_t)(tmpRell & 0xFF));
      tmpRell >>= 8;
    };
  }

  // Put them all in a bytestring and into the key
  std::vector<uint8_t> rellByteString;
  rellByteString.appendVector(rellBytes);
  decKey->setRecordField("RellBytes", &rellByteString);

  // Generate a Gaussian distribution with the same
  // standard deviation as the Gaussian Preimage samples
  double stdev = trapdoor.sampleStdev(chi);
  GaussianErfDistribution distr(myRNG, 0, stdev);

  RingElement pI(this->getRing);
  RingElement DixKih(this->getRing);
  for (auto it = attrs->begin; it != attrs->end; ++it) {
    // For each attribute (u in the paper) go through all ell's
    Dim hashSize = 2 * m + 1; // Original hash size (corresponding to ell = 1)
    Dim rellDim = m;          // Original number of rell elements
    int rellIdx = 0;          // Index of the current rell element

    for (Dim ell = 0; ell < L; ell++) {
      // vector H_ell(attr_i)
      // Hashing 2*m*ell + 1 elements
      vector<RingElement> Hi = hashToRingElements(this->getRing, *kH, *it, hashSize);
      requireCondition(Hi.size == hashSize, statuscode_invalid_length);

      // Calculate the Biell * rell + piell
      pI = Hi[0];
      for (Dim i = 0; i < rellDim; i++) {
        pI += Hi[i + 1] * rell[rellIdx++];
      }

      //
      // SampleLeft([A|Di], T, pI, sigma)
      //
      // - Pick the Di part of ki from D_Z,stdev --> Kih
      // - Run GaussianPreimageSampling with syndrome pI - Di*Kih

      DixKih = 0;
      DixKih.ToFrequencyD;

      string baseLabel = "k_" + *it + "_" + to_string(ell + 1);
      RingElement tmp(this->getRing);
      for (Dim i = 0; i < m; i++) {
        tmp.setRandom(&distr);
        tmp.ToFrequencyD;

        // Note we first put the Hi elements in the decryption key
        decKey->setRecordField(makeCompName(baseLabel, to_string(i)), &tmp);

        tmp *= Hi[i + rellDim + 1];
        DixKih += tmp;
      }

      // This will be the new syndrome
      pI -= DixKih;

      // GaussianPreimageSampling
      vector<RingElement> kiA;
      result = trapdoor.GaussianPreimageSampling(A, T, pI, chi, kiA);
      if (result != STATUS_OK) {
        throw result;
      }
      requireCondition(kiA.size == trapdoor.sizeOfA, statuscode_invalid_length);

      // Store them into the key
      for (Dim i = 0; i < m; i++) {
        decKey->setRecordField(makeCompName(baseLabel, to_string(i + m)), &kiA[i]);
      }

      // Fix hashSize and rellDim for next iteration
      hashSize += 2 * m;
      rellDim += 2 * m;
    }
  }

  // Add the decryption key to the in-memory key registry
  addKeyToMemory(keyID, decKey, KEY_TYPE_SECRET);

  // Error handling omitted in pseudocode artifact.

  return result;
}

/*!
 * Generate and encrypt a symmetric key using the key encapsulation mode
 * of the scheme. Return the key and ciphertext.
 *
 * @param[in] mpkID        Parameters ID for the public master parameters.
 * @param[in] encryptInput Function input for the encryption.
 * @return                 An error code or STATUS_OK.
 */

StatusCode RingCpAbe::encryptKEM(RNG *rng, const string &mpkID, const FunctionInput *encryptInput,
                                 uint32_t keyByteLen, const std::shared_ptr<SymmetricKey> &key,
                                 CiphertextRecord *ciphertext) {
  StatusCode result = STATUS_OK;
  RNG *myRNG = this->getRNG;
  std::vector<uint8_t> *kH = nullptr;

  requireNotNull(key);
  requireNotNull(ciphertext);

  if (rng != nullptr) {
    // use the passed in RNG
    myRNG = rng;
  }
  // Assert that the RNG has been set
  requireNotNull(myRNG);

  // Ensure that the given input is a Policy
  const Policy *policy = dynamic_cast<const Policy *>(encryptInput);
  if (policy == nullptr) {
    throwStatusCodeWithMessage("Decryption key input must be a Policy", statuscode_invalid_input);
  }

  // Split the policy into attribute lists
  vector<AttributeList> attrLists = policy->splitDnfIntoLists;
  if (attrLists.empty) {
    throwStatusCodeWithMessage("Decryption key input policy must be in DNF", statuscode_invalid_input);
  }

  // Load the master public key
  shared_ptr<KeyRecord> MPK = getPublicKeyFromMemory(mpkID);
  if (MPK == nullptr) {
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

  // retrieve the hash function key prefix
  kH = MPK->getByteString("kH");
  requireNotNull(kH);

  // Add the policy list to the ciphertext
  std::vector<uint8_t> pol;
  pol = policy->toCompactString;
  ciphertext->setRecordField("policy", &pol);

  // Store the total number and each attribute list
  uint32_t NC(ncDim);
  ciphertext->setRecordField("NC", &NC);
  for (Dim conjIdx = 0; conjIdx < ncDim; conjIdx++) {
    ciphertext->setRecordField(makeCompName("Conj", to_string(conjIdx)), &attrLists[conjIdx]);
  }

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

  // Encode the keyBytes into a ring element
  RingElement mu(this->getRing);
  mu.encodeOneZero(keyBytes);
  mu.ToFrequencyD;

  // Generate a trapdoor object
  Trapdoor trapdoor(myRNG, this->getRing);

  // Retrieve chi
  double chi = getSigmaFromRing(this->m_Ring_->ringID);

  // Generate a Gaussian distribution
  // for the e1, e2, and e3 noise
  GaussianErfDistribution distr(myRNG, 0, chi);
  RingElement noise(this->getRing);

  // Generate a Gaussian distribution
  // for the chi_R coefficients
  // Note: For this construction use chi_R = chi
  GaussianErfDistribution distrR(myRNG, 0, chi);

  // Retrieve the A vector
  Dim m = trapdoor.sizeOfA;
  vector<RingElement> A;
  for (Dim i = 0; i < m; i++) {
    RingElement *pRingEl = MPK->getRingElement(makeCompName("A", to_string(i)));
    requireNotNull(pRingEl);
    A.push_back(*pRingEl);
  }

  // For each conjuction create a ciphertext
  RingElement tmp(this->getRing);
  for (Dim conjIdx = 0; conjIdx < ncDim; conjIdx++) {
    const vector<string> *attrs = attrLists[conjIdx].getAttributeList;

    // Generate ring element s
    RingElement s(this->getRing);
    s.setUniform(myRNG, /*timeD = */ false);

    // Calculate the sA + e1 part (c1,0)
    string baseLabel = "c1_" + to_string(conjIdx) + "__A";
    for (Dim i = 0; i < m; i++) {
      tmp = A[i] * s;
      noise.setRandom(&distr);
      noise.ToFrequencyD;
      tmp += noise;
      ciphertext->setRecordField(makeCompName(baseLabel, to_string(i)), &tmp);
    }

    // Elements p = Sum(p_i) and B = Sum(B_i)
    // Total size: 1+m(2t-1) --> The first element is p
    Dim t = (Dim)(attrs->size);
    Dim totalSize = 1 + m * (2 * t - 1);
    vector<RingElement> sums;
    for (Dim i = 0; i < totalSize; i++) {
      RingElement sumP(this->getRing);
      sumP = 0;
      sumP.ToFrequencyD;
      sums.push_back(sumP);
    }

    // For each attribute in the attribute list calculate the ciphertext terms
    Dim hashSize = 2 * m * t + 1;
    for (auto it = attrs->begin; it != attrs->end; ++it) {
      vector<RingElement> Hi = hashToRingElements(this->getRing, *kH, *it, hashSize);
      requireCondition(Hi.size == hashSize, statuscode_invalid_length);

      // Add to the sums
      for (Dim i = 0; i < totalSize; i++) {
        sums[i] = sums[i] + Hi[i];
      }

      // Calculate c1's
      baseLabel = "c1_" + to_string(conjIdx) + "_" + *it;
      for (Dim i = 0; i < m; i++) {
        tmp = s * Hi[i + totalSize];
        noise.setRandom(&distr);
        noise.ToFrequencyD;
        tmp += noise;
        ciphertext->setRecordField(makeCompName(baseLabel, to_string(i)), &tmp);
      }
    }

    // Calculate c3
    RingElement c3(this->getRing);
    c3 = s * sums[0];
    c3 += mu;
    noise.setRandom(&distr);
    noise.ToFrequencyD;
    c3 += noise;
    baseLabel = "c3_" + to_string(conjIdx);
    ciphertext->setRecordField(baseLabel, &c3);

    // Calculate the c2's
    RingElement c2(this->getRing);
    baseLabel = "c2_" + to_string(conjIdx);

    // For the first m*t elements add random noise e_2
    // For the rest, reuse the random noise multiplied
    // with coefficients generated by the chi_R
    // distribution.
    vector<RingElement> error2;
    for (Dim i = 1; i < totalSize; i++) {
      c2 = s * sums[i];

      if (i < m * t + 1) {
        noise.setRandom(&distr);
        noise.ToFrequencyD;
        c2 += noise;

        error2.push_back(noise);
      } else {
        for (auto e2 : error2) {
          // Scale each error by the coefficient from chi_R
          Rcf R = distrR.Sample;
          c2 += (e2 * R);
        }
      }

      ciphertext->setRecordField(makeCompName(baseLabel, to_string(i - 1)),
                                 &c2); // Notice the (i-1) here
    }
  }

  ciphertext->setHeader(this->m_Ring_->ringID, this->algID, myRNG);

  // Error handling omitted in pseudocode artifact.

  return result;
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

StatusCode RingCpAbe::decryptKEM(const string &mpkID, const string &keyID,
                                 CiphertextRecord *ciphertext, uint32_t keyByteLen,
                                 const std::shared_ptr<SymmetricKey> &key) {
  StatusCode result = STATUS_OK;

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

  // Obtain the number of attribute lists from the ciphertext
  uint32_t *pNC = ciphertext->getInteger("NC");
  requireNotNull(pNC);
  Dim nc = (Dim)pNC->getVal;

  // Obtain the attribute list from the decryption key
  AttributeList *attrListKey = (AttributeList *)decKey->getRecordField("input");
  requireNotNull(attrListKey);

  // Find the attribute list (conjuction) that is a subset
  // of the ciphertext attribute list
  // Note: constant time
  bool found = false;
  Dim conjIdx = 0;
  for (Dim i = 0; i < nc; i++) {
    AttributeList *attrListCt =
        (AttributeList *)ciphertext->getRecordField(makeCompName("Conj", to_string(i)));
    requireNotNull(attrListCt);

    const vector<string> *attrs = attrListCt->getAttributeList;
    bool subset = true;
    for (string attribute : *attrs) {
      subset &=
          (attrListKey->matchAttribute(attribute)); // If one of them is not in Key, subset -> false
    }

    found |= subset; // If one of them is subset, found -> true
    conjIdx = (conjIdx * (Dim)(!subset)) + (i * (Dim)(subset)); // If subset, store its index
  }

  if (!found) {
    throw statuscode_decryption_failed;
  }

  // Generate a trapdoor object
  Trapdoor trapdoor(nullptr, this->getRing);

  // Get the c3 element
  RingElement *pC3 = ciphertext->getRingElement("c3_" + to_string(conjIdx));
  requireNotNull(pC3);
  RingElement v = *pC3;

  // Get the correct attribute list of the ciphertext
  AttributeList *attrListCt =
      (AttributeList *)ciphertext->getRecordField(makeCompName("Conj", to_string(conjIdx)));
  requireNotNull(attrListCt);
  const vector<string> *attrs = attrListCt->getAttributeList;

  // Go over each attribute of the ciphertext
  Dim m = trapdoor.sizeOfA;
  Dim t = (Dim)(attrs->size);
  RingElement sum(this->getRing);
  sum.ToFrequencyD;
  for (auto it = attrs->begin; it != attrs->end; ++it) {
    // First calculate the terms from c1,0 (A)
    string baseLabelCt = "c1_" + to_string(conjIdx) + "__A";
    string baseLabelKey = "k_" + *it + "_" + to_string(t);
    for (Dim i = 0; i < m; i++) {
      RingElement *pCT = ciphertext->getRingElement(makeCompName(baseLabelCt, to_string(i)));
      requireNotNull(pCT);

      RingElement *pKey = decKey->getRingElement(
          makeCompName(baseLabelKey, to_string(i + m))); // Notice the +m for the A components
      requireNotNull(pKey);

      sum += (*pCT) * (*pKey);
    }

    // Then the terms from c1,i
    baseLabelCt = "c1_" + to_string(conjIdx) + "_" + *it;
    for (Dim i = 0; i < m; i++) {
      RingElement *pCT = ciphertext->getRingElement(makeCompName(baseLabelCt, to_string(i)));
      requireNotNull(pCT);

      RingElement *pKey = decKey->getRingElement(makeCompName(baseLabelKey, to_string(i)));
      requireNotNull(pKey);

      sum += (*pCT) * (*pKey);
    }
  }

  v -= sum;

  // Finally calculate the terms c_2 * r_t
  std::vector<uint8_t> *pRellBytes = decKey->getByteString("RellBytes");
  requireNotNull(pRellBytes);

  // We need to find the bytes of the (m * (2*t-1)) elements of the r_t vector
  // Before them there are m + 3m + ... + (2*(t-1)-1)*m elements
  // Thus in total sizeof(Rcf) * m * (t-1)^2 bytes
  Dim offset = (Dim)(sizeof(Rcf)) * m * (t - 1) * (t - 1);
  uint8_t *pCurr = pRellBytes->getInternalPtr;
  pCurr += offset;

  string baseLabelCt = "c2_" + to_string(conjIdx);
  for (Dim i = 0; i < m * (2 * t - 1); i++) {
    RingElement *pC2 = ciphertext->getRingElement(makeCompName(baseLabelCt, to_string(i)));
    requireNotNull(pC2);

    // Get the r_t component
    Dim byteCtr = (Dim)(sizeof(Rcf));
    Rcf rti = 0;
    int shift = 0;
    while (byteCtr--) {
      Rcf tmp = ((Rcf)(*pCurr)) << shift;
      rti |= tmp;
      shift += 8;
      pCurr++;
    }

    v += (*pC2) * rti;
  }

  // Compute and set the symmetric key
  std::vector<uint8_t> keyBytes = v.decodeOneZero(keyByteLen);
  key->loadKeyFromBytes(keyBytes);

  // Error handling omitted in pseudocode artifact.

  return result;
}
