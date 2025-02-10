#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (rows <= 0 || columns <= 0 || result == NULL) {
    return INCORRECT;
  }
  int result_code = OK;
  result->rows = rows;
  result->columns = columns;
  result->matrix = calloc(rows, sizeof(double *));
  if (!result->matrix) {
    result_code = INCORRECT;
  } else {
    for (int i = 0; i < rows && result_code == OK; ++i) {
      result->matrix[i] = calloc(columns, sizeof(double));
      if (!result->matrix[i]) {
        while (i--) {
          free(result->matrix[i]);
        }
        free(result->matrix);
        result->matrix = NULL;
        result_code = INCORRECT;
      }
    }
  }
  return result_code;
}

void s21_remove_matrix(matrix_t *A) {
  if (A && A->matrix) {
    for (int i = 0; i < A->rows; ++i) {
      if (A->matrix[i]) {
        free(A->matrix[i]);
        A->matrix[i] = NULL;
      }
    }
    free(A->matrix);
    A->matrix = NULL;
  }
  if (A) {
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (A->rows != B->rows || A->columns != B->columns || A->rows <= 0 ||
      B->rows <= 0) {
    return FAILURE;
  }
  int compare_result = SUCCESS;
  for (int i = 0; i < A->rows && compare_result; ++i) {
    for (int j = 0; j < A->columns && compare_result; ++j) {
      if (!s21_compare_double(A->matrix[i][j], B->matrix[i][j])) {
        compare_result = FAILURE;
      }
    }
  }
  return compare_result;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int result_code = OK;
  s21_arguments_check(A, B, &result_code, NOT_MUL_MATRIX);
  if (result_code == OK) {
    result_code = s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows && result_code == OK; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }
  return result_code;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int result_code = OK;
  s21_arguments_check(A, B, &result_code, NOT_MUL_MATRIX);
  if (result_code == OK) {
    result_code = s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows && result_code == OK; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return result_code;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int result_code = OK;
  s21_arguments_check(A, A, &result_code, NOT_MUL_MATRIX);
  if (result_code == OK) {
    if (isinf(number) || isnan(number)) {
      result_code = CALC_ERROR;
    } else {
      result_code = s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows && result_code == OK; ++i) {
        for (int j = 0; j < A->columns; ++j) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    }
  }
  return result_code;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int result_code = OK;
  s21_arguments_check(A, B, &result_code, MUL_MATRIX);
  if (result_code == OK) {
    result_code = s21_create_matrix(A->rows, B->columns, result);
    double *A_row = calloc(A->columns, sizeof(double));
    double *B_column = calloc(B->rows, sizeof(double));
    if (A_row && B_column) {
      for (int row = 0; row < A->rows && result_code == OK; ++row) {
        s21_row_copy(A, A_row, row);
        for (int column = 0; column < B->columns; ++column) {
          s21_column_copy(B, B_column, column);
          result->matrix[row][column] = s21_row_column_mul(A_row, B_column, A);
        }
      }
    } else {
      result_code = INCORRECT;
    }
    free_line(A_row);
    free_line(B_column);
  }
  return result_code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int result_code = OK;
  s21_arguments_check(A, A, &result_code, NOT_MUL_MATRIX);
  if (result_code == OK) {
    result_code = s21_create_matrix(A->columns, A->rows, result);
    if (result_code == OK) {
      s21_transpose_additional(A, result);
    }
  }
  return result_code;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int result_code = OK;
  s21_arguments_check(A, A, &result_code, DETERMINANT);
  if (result_code == OK) {
    result_code = s21_create_matrix(A->rows, A->columns, result);
    if (result_code == OK) {
      result_code = s21_calc_complements_additional(A, result);
    }
  }
  return result_code;
}

int s21_determinant(matrix_t *A, double *result) {
  int result_code = OK;
  s21_arguments_check(A, A, &result_code, DETERMINANT);
  if (result_code == OK) {
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else {
      *result = s21_calc_determinant(A, A->rows);
    }
  }
  return result_code;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int result_code = OK;
  double determinant = 0.0, inverse_det = 0.0;
  s21_arguments_check(A, A, &result_code, INVERSE);
  if (result_code == OK) {
    determinant = s21_calc_determinant(A, A->rows);
    if (s21_compare_double(determinant, 0.0)) {
      result_code = CALC_ERROR;
    } else {
      result_code = s21_create_matrix(A->rows, A->columns, result);
      if (result_code == OK) {
        matrix_t temp = {0};
        if (s21_create_matrix(A->rows, A->columns, &temp) == OK) {
          inverse_det = 1.0 / determinant;
          s21_calc_complements_additional(A, &temp);
          s21_transpose_additional(&temp, result);
          s21_remove_matrix(&temp);
          for (int i = 0; i < result->rows; ++i) {
            for (int j = 0; j < result->columns; ++j) {
              result->matrix[i][j] *= inverse_det;
            }
          }
        }
      }
    }
  }
  return result_code;
}
