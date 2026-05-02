#include "ring.h"

///
/// \file   ring_presets.cpp
///
/// \brief  Implementation of the PresetRings out-of-class definition.
///

#include "lwe_base_types.h"
#include "ring.h"
#include "ring_element.h"
#include "field_element.h"
#include "distribution.h"
#include "trapdoor.h"
#include "kpabe.h"
#include "cpabe.h"


    // This is the out-of-class definition. The parameters are initialized in ring.h
    constexpr std::array<RingParameters, 6> PresetRings::Parameters;

#ifndef __RCF64__
    std::ostream& operator<<(std::ostream& os, Rcf value) {
        URcf tmp = value < 0 ? -value : value;
        char buffer[129];           // It fits the entire value + the null terminated char
        char* d = std::end(buffer);

        // Put the null terminator
        --d;
        *d = '\0';

        do {
            -- d;
            *d = "0123456789"[tmp % 10];
            tmp /= 10;
        } while (tmp != 0);

        if (value < 0) {
            -- d;
            *d = '-';
        }

        std::string number(d);
        os << number;

        return os;
    }
#endif

    size_t decimalLength(Rcf value) {
        size_t res = 0;

        if (value < 0) {
            res++;
        }

        do {
            res++;
            value /= 10;
        } while (value != 0);

        return res;
    }

