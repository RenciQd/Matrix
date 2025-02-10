#ifndef TEST_H
#define TEST_H

#include <time.h>
#include <unistd.h>

#include "../s21_matrix.h"
#include "check.h"

Suite *suite_create_matrix(void);
Suite *suite_eq_matrix(void);
Suite *suite_sum_matrix(void);
Suite *suite_sub_matrix(void);
Suite *suite_mult_number_matrix(void);
Suite *suite_mult_matrix(void);
Suite *suite_transpose_matrix(void);
Suite *suite_calc_complements(void);
Suite *suite_determinant(void);
Suite *suite_inverse_matrix(void);
Suite *suite_matrix_gather(void);

void run_testcase(Suite *testcase);
double get_rand(double min, double max);
void run_tests(void);

#endif