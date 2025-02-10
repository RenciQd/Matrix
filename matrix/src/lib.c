#include "s21_matrix.h"

int s21_compare_double(double first, double second) {
  return (fabs(first - second) > EPSILON) ? NOT_EQUAL : EQUAL;
}

int s21_nan_inf_check(matrix_t *value) {
  int has_invalid_value = 0;
  for (int i = 0; i < value->rows && !has_invalid_value; i++) {
    for (int j = 0; j < value->columns && !has_invalid_value; j++) {
      if (isnan(value->matrix[i][j]) || isinf(value->matrix[i][j])) {
        has_invalid_value = 1;
      }
    }
  }
  return has_invalid_value;
}

void s21_arguments_check(matrix_t *value_1, matrix_t *value_2, int *result_code,
                         int type) {
  if (s21_nan_inf_check(value_1) || s21_nan_inf_check(value_2)) {
    *result_code = INCORRECT;
  } else if (!s21_is_size_ok(value_1, value_2) || value_1 == NULL ||
             value_2 == NULL) {
    *result_code = INCORRECT;
  }
  if (*result_code == OK) {
    if (type == NOT_MUL_MATRIX && (value_1->columns != value_2->columns ||
                                   value_1->rows != value_2->rows)) {
      *result_code = CALC_ERROR;
    } else if (type == MUL_MATRIX && (value_1->rows != value_2->columns ||
                                      value_1->columns != value_2->rows)) {
      *result_code = CALC_ERROR;
    } else if ((type == DETERMINANT || type == INVERSE) &&
               value_1->rows != value_1->columns) {
      *result_code = CALC_ERROR;
    }
  }
}

int s21_is_size_ok(matrix_t *value_1, matrix_t *value_2) {
  return (value_1->rows > 0 && value_1->columns > 0 && value_2->rows > 0 &&
          value_2->columns > 0)
             ? SIZE_OK
             : SIZE_NOT_OK;
}

void s21_matrix_copy(matrix_t *source, matrix_t *destination) {
  for (int i = 0; i < source->rows; i++) {
    for (int j = 0; j < source->columns; j++) {
      destination->matrix[i][j] = source->matrix[i][j];
    }
  }
}

void s21_row_copy(matrix_t *matrix, double *row_data, int row_index) {
  for (int i = 0; i < matrix->columns; i++) {
    row_data[i] = matrix->matrix[row_index][i];
  }
}

void s21_column_copy(matrix_t *matrix, double *column_data, int column_index) {
  for (int i = 0; i < matrix->rows; i++) {
    column_data[i] = matrix->matrix[i][column_index];
  }
}

double s21_row_column_mul(double *row_data, double *column_data,
                          matrix_t *matrix) {
  double result = 0.0;
  for (int i = 0; i < matrix->columns; i++) {
    result += row_data[i] * column_data[i];
  }
  return result;
}

void s21_get_cofactor(matrix_t *source, matrix_t *minor, int excluded_row,
                      int excluded_column, int dimension) {
  int row_pos = 0, col_pos = 0;
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      if (i != excluded_row && j != excluded_column) {
        minor->matrix[row_pos][col_pos++] = source->matrix[i][j];
        if (col_pos == minor->rows) {
          col_pos = 0;
          row_pos++;
        }
      }
    }
  }
}

double s21_calc_determinant(matrix_t *matrix, int size) {
  double sign = 1.0, determinant = 0.0;
  matrix_t temp_minor = {0};
  if (size == 1) return matrix->matrix[0][0];
  s21_create_matrix(size - 1, size - 1, &temp_minor);
  for (int col = 0; col < size; col++) {
    s21_get_cofactor(matrix, &temp_minor, 0, col, size);
    determinant += sign * matrix->matrix[0][col] *
                   s21_calc_determinant(&temp_minor, size - 1);
    sign = -sign;
  }
  s21_remove_matrix(&temp_minor);
  return determinant;
}

void s21_transpose_additional(matrix_t *matrix, matrix_t *transposed) {
  for (int i = 0; i < matrix->rows; i++) {
    for (int j = 0; j < matrix->columns; j++) {
      transposed->matrix[j][i] = matrix->matrix[i][j];
    }
  }
}

int s21_calc_complements_additional(matrix_t *matrix, matrix_t *complements) {
  int status = OK;
  if (matrix->rows == 1) {
    complements->matrix[0][0] = 1;
  } else {
    double minor_value = 0.0;
    matrix_t minor_matrix = {0};
    status =
        s21_create_matrix(matrix->rows - 1, matrix->columns - 1, &minor_matrix);
    for (int i = 0; i < matrix->rows && matrix->rows > 1 && status == OK; i++) {
      for (int j = 0; j < matrix->columns; j++) {
        s21_get_cofactor(matrix, &minor_matrix, i, j, matrix->rows);
        minor_value = s21_calc_determinant(&minor_matrix, minor_matrix.columns);
        complements->matrix[i][j] = minor_value;
      }
    }
    s21_negate_matrix(complements);
    s21_remove_matrix(&minor_matrix);
  }
  return status;
}

void free_line(double *line_data) {
  if (line_data != NULL) {
    free(line_data);
    line_data = NULL;
  }
}

void s21_negate_matrix(matrix_t *matrix) {
  for (int i = 0; i < matrix->rows; i++) {
    for (int j = 0; j < matrix->columns; j++) {
      matrix->matrix[i][j] *= ((i + j) % 2) ? -1.0 : 1.0;
    }
  }
}