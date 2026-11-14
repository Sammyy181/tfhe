cmake_minimum_required(VERSION 3.0)

# Toom-4 striding transform processor (Bodrato interpolation), exact replacement for the FFT

set(SRCS
    lagrangehalfc_impl.cpp
    )

set(HEADERS
    lagrangehalfc_impl.h
    )

if (ENABLE_TOOM4_STRIDED)
    add_library(tfhe-fft-toom4-strided OBJECT ${SRCS} ${HEADERS})
    set_property(TARGET tfhe-fft-toom4-strided PROPERTY POSITION_INDEPENDENT_CODE ON)
endif (ENABLE_TOOM4_STRIDED)
