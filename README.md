# Large-Universe (Multi-Authority) ABE from LWE

Pseudocode artifact for CCS 2026 (Cycle B), covering Ring-LWE KP-ABE and CP-ABE.

## Scope
- In scope: Ring-LWE KP-ABE and CP-ABE .
- Out of scope: integer-LWE, multi-authority variants, short-key variants.

## What this is
Library-neutral pseudocode adaptation for reviewer inspection. It omits generic components for parsing, satisfiability checks, secret sharing, hybrid encryption, and key serialization/storage. The code (with the exception of `tools/prime_moduli.cpp`) does not compile.

## Files
- `include/*`: base types, ring, elements, distributions, trapdoor, scheme APIs.
- `src/*`: pseudocode implementations for ring utilities and KP/CP algorithms.
- `tools/prime_moduli.cpp`: standalone parameter-generation tool.

## Mapping to paper
- `ring.*`: ring construction and presets.
- `distribution.*`: RNG and distribution helpers.
- `trapdoor.*`: trapdoor generation and sampling.
- `kpabe.*`: KP-ABE setup/keygen/encrypt/decrypt pseudocode.
- `cpabe.*`: CP-ABE setup/keygen/encrypt/decrypt pseudocode.
