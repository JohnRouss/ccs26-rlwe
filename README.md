# Large-Universe (Multi-Authority) ABE from LWE
Source code for the RingLWE constructions for the CCS 2026 (Cycle B) submission.

## Header Files
Directory: include

| File	| Description	| Scheme|
|:--------|:----------|:----------|
| contextcp.h	| Ring-LWE CP-ABE scheme API	| CP-ABE|
| contextkp.h	| Ring-LWE KP-ABE scheme API	| KP-ABE|
| distribution.h	| Header file for the different distributions used by the PQC schemes (Uniform and Gaussians)| 	Both|
| fieldelement.h	| Field element class, operations on field elements and their spectra (APIs)| 	Both|
| lwebasetypes.h	| Base types for the ring coefficients and LWE constants. Also, implementations of the low-level mathematical operations.| 	Both|
| ring.h	| Ring class and parameters| 	Both|
| ringelement.h	| Ring element class and operations on ring elements| 	Both|
| trapdoor.h	| Trapdoor generation API and pre-image sampling API| 	Both|

## Implementation files
Directory: src
| File | Description | Scheme |
|:--------|:----------|:----------|
| contextcp.cpp	| Ring-LWE CP-ABE scheme	| CP-ABE| 
| contextkp.cpp	| Ring-LWE KP-ABE scheme	| KP-ABE| 
| distribution.cpp	| Implementation for the different distributions used by the PQC schemes (Uniform and Gaussians)| Both| 
| fieldelement.cpp	| Field element and spectra operations| Both| 
| ring.cpp	| Constructors and setup algorithms for rings| Both| 
| ringelement.cpp	| Ring element operations| Both| 
| ringpresets.cpp	| Some helper objects for ring parameters and printing ring elements| Both| 
| trapdoor.cpp	| Trapdoor generation and pre-image sampling implementations| Both| 

## Tools
Directory: tools
| File | Description | Scheme |
|:--------|:----------|:----------|
| prime_moduli.cpp	| Script file to generate and print the prime moduli and parameters needed for the Ring LWE constructions.	| Both| 
