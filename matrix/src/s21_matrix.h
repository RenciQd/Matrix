#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define OK 0
#define INCORRECT 1  // Ошибка, некорректная матрица;

// Ошибка вычисления (несовпадающие размеры матриц; матрица,
// для которой нельзя провести вычисления и т. д.).
#define CALC_ERROR 2
#define SUCCESS 1
#define FAILURE 0
#define EQUAL 1
#define NOT_EQUAL 0

#define EPSILON 9e-8
// macros for s21_arguments_check
#define NOT_MUL_MATRIX 1
#define MUL_MATRIX 2
#define DETERMINANT 3
#define INVERSE 4
#define COMPLEMENT 5
// for s21_is_size_ok
#define SIZE_NOT_OK 0
#define SIZE_OK 1

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

// matrix.c
int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

// lib.c
int s21_compare_double(double first, double second);
int s21_nan_inf_check(matrix_t *value);
// void s21_nan_handling(int *result_code, matrix_t *value);
void s21_arguments_check(matrix_t *value_1, matrix_t *value_2, int *result_code,
                         int type);
int s21_is_size_ok(matrix_t *value_1, matrix_t *value_2);
void s21_matrix_copy(matrix_t *value, matrix_t *result);
void s21_row_copy(matrix_t *value_1, double *A_row, int row);
void s21_column_copy(matrix_t *B, double *B_column, int column);
double s21_row_column_mul(double *A_row, double *B_column, matrix_t *A);
void s21_get_cofactor(matrix_t *value, matrix_t *temp, int row, int column,
                      int size);
double s21_calc_determinant(matrix_t *value, int size);
void s21_transpose_additional(matrix_t *value, matrix_t *result);
int s21_calc_complements_additional(matrix_t *value, matrix_t *result);
void free_line(double *value);
void s21_negate_matrix(matrix_t *value);

// helpers.c
void print_line(double *line, int elements);
void print_matrix(matrix_t matrix);