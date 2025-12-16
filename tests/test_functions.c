/**
 * @file test_functions.c
 * @brief Unit tests for mathematical and utility functions
 * @details This test suite validates the functionality of various mathematical operations
 *          and random number generation functions used in the TFM project.
 *
 * @test test_sqrt_distance()
 *       Verifies the Euclidean distance calculation for various input cases:
 *       - Standard case (3,4) should equal 5
 *       - Zero vector should return 0
 *       - Unit diagonal (1,1) should equal sqrt(2)
 *
 * @test test_scalar_product()
 *       Validates dot product computation:
 *       - (1,2)·(3,4) should equal 11
 *       - Zero vector cases should return 0
 *       - Orthogonal vectors should return 0
 *
 * @test test_PBC()
 *       Tests periodic boundary conditions with unit cell:
 *       - Positive values > 0.5 wrap to negative range
 *       - Negative values < -0.5 wrap to positive range
 *       - Zero remains unchanged
 *
 * @test test_inverse_cumulative_rayleigh()
 *       Validates inverse cumulative distribution function for Rayleigh distribution
 *       - Result must be positive and bounded
 *
 * @test test_new_vector_segment()
 *       Verifies vector component generation from magnitude and angle:
 *       - Angle 0 produces full x-component, zero y-component
 *       - Angle π/2 produces zero x-component, full y-component
 *
 * @test test_calculate_f()
 *       Tests scalar field evaluation at specified indices
 *
 * @test test_random_functions()
 *       Validates random number generation within specified ranges:
 *       - ini_ran() initializes the RNG with seed
 *       - randomInPR() generates values within given bounds
 *
 * @test test_box_muller()
 *       Verifies Gaussian random number generation produces finite results
 *
 * @dependencies head.h, math.h, stdlib.h, stdio.h, assert.h
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "head.h"

// Test sqrt_distance
void test_sqrt_distance() {
    assert(fabs(sqrt_distance(3.0, 4.0) - 5.0) < 1e-9);
    assert(fabs(sqrt_distance(0.0, 0.0) - 0.0) < 1e-9);
    assert(fabs(sqrt_distance(1.0, 1.0) - sqrt(2.0)) < 1e-9);
    printf("OK sqrt_distance tests passed\n");
}

// Test scalar_product
void test_scalar_product() {
    assert(fabs(scalar_product(1.0, 2.0, 3.0, 4.0) - 11.0) < 1e-9);
    assert(fabs(scalar_product(0.0, 0.0, 5.0, 5.0) - 0.0) < 1e-9);
    assert(fabs(scalar_product(1.0, 0.0, 0.0, 1.0) - 0.0) < 1e-9);
    printf("OK scalar_product tests passed\n");
}

// Test PBC
void test_PBC() {
    double x = 0.6;
    PBC(&x, 1.0);
    assert(fabs(x - (-0.4)) < 1e-9);

    x = -0.6;
    PBC(&x, 1.0);
    assert(fabs(x - (0.4)) < 1e-9);

    x = 0.0;
    PBC(&x, 1.0);
    assert(fabs(x - 0.0) < 1e-9);

    printf("OK PBC tests passed\n");
}

// Test inverse_cumulative_rayleigh
void test_inverse_cumulative_rayleigh() {
    double result = inverse_cumulative_rayleigh(0.5, 1.0);
    assert(result > 0.0);
    assert(result < 2.0);
    printf("OK inverse_cumulative_rayleigh tests passed\n");
}

// Test new_vector_segment
void test_new_vector_segment() {
    double dx, dy;
    new_vector_segment(1.0, 0.0, &dx, &dy);
    assert(fabs(dx - 1.0) < 1e-9);
    assert(fabs(dy - 0.0) < 1e-9);
    
    new_vector_segment(1.0, PI / 2.0, &dx, &dy);
    assert(fabs(dx - 0.0) < 1e-9);
    assert(fabs(dy - 1.0) < 1e-9);
    printf("OK new_vector_segment tests passed\n");
}

// Test calculate_f
void test_calculate_f() {
    double x[] = {0.0, 1.0};
    double y[] = {0.0, 1.0};
    double f = calculate_f(x, y, 0, 1.0, 1.0);
    assert(fabs(f - 0.0) < 1e-9);
    
    f = calculate_f(x, y, 1, 1.0, 1.0);
    assert(fabs(f - 1.0) < 1e-9);
    printf("OK calculate_f tests passed\n");
}

// Test ini_ran and randomInPR
void test_random_functions() {
    ini_ran(12345);
    float r1 = randomInPR(0.0, 1.0);
    assert(r1 >= 0.0 && r1 <= 1.0);
    
    float r2 = randomInPR(5.0, 10.0);
    assert(r2 >= 5.0 && r2 <= 10.0);
    printf("OK random functions tests passed\n");
}

// Test box_muller
void test_box_muller() {
    ini_ran(54321);
    double result = box_muller();
    assert(isfinite(result));
    printf("OK box_muller tests passed\n");
}

int main() {
    test_sqrt_distance();
    test_scalar_product();
    test_PBC();
    test_inverse_cumulative_rayleigh();
    test_new_vector_segment();
    test_calculate_f();
    test_random_functions();
    test_box_muller();
    
    printf("\nOK All tests passed!\n");
    return 0;
}