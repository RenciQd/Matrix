CHECK_CFLAGS := -I/opt/homebrew/opt/check/include
CHECK_LIBS := -L/opt/homebrew/opt/check/lib -lcheck
FLAGS = -D_POSIX_C_SOURCE=200809L -Wall -Werror -Wextra -std=c11
GCOV_FLAGS = -fprofile-arcs -ftest-coverage

TEST_FILES = ./tests/*.c
C_FILES = lib.c s21_matrix.c
H_FILES = s21_matrix.h

OS := $(shell uname -s)
ISWSL := $(shell uname -r | grep WSL -o)

ifeq ($(shell uname -r | grep WSL -o), WSL)
	HTML_OPEN = @echo HTML report is in $(shell pwd)/report/index.html
else ifeq ($(shell uname -s), Linux)
	HTML_OPEN = @xdg-open ./report/index.html
else
	HTML_OPEN = @open ./report/index.html
endif

all: s21_matrix.a test gcov_report

s21_matrix.a:
	gcc $(FLAGS) -c -o s21_matrix.o s21_matrix.c
	gcc $(FLAGS) -c -o lib.o lib.c
	ar rcs s21_matrix.a s21_matrix.o lib.o
	@rm -f lib.o s21_matrix.o

gcov_report:
	@rm -rf report
	gcc $(FLAGS) $(GCOV_FLAGS) $(C_FILES) $(TEST_FILES) -o gcov_report.out $(CHECK_CFLAGS) $(CHECK_LIBS)
	./gcov_report.out
	lcov -t "s21_matrix" -o s21_matrix.info -c -d . -q
	genhtml -o report s21_matrix.info -q
	@rm -f s21_matrix.info *.gcda *.gcno gcov_report.out
	$(HTML_OPEN)

test: s21_matrix.a
	gcc -o test_out $(FLAGS) $(TEST_FILES) s21_matrix.a $(CHECK_CFLAGS) $(CHECK_LIBS)
	./test_out
	rm *.a *.out
	rm -rf report

rebuild: clean all

clean:
	rm -rf report
	rm -f *.a test_out *.out *.o *.h.gch *.info

clang:
	clang-format -i --style=Google $(TEST_FILES) $(C_FILES) $(H_FILES)

memory_test_gcc_s21_matrix.a:
	gcc $(FLAGS) -c -o s21_matrix.o s21_matrix.c -lm  -fsanitize=address
	gcc $(FLAGS) -c -o lib.o lib.c -lm -fsanitize=address
	ar rcs s21_matrix.a s21_matrix.o lib.o
	@rm -f s21_matrix.o lib.o

memory_test_tests: s21_matrix.a
	gcc $(FLAGS) $(TEST_FILES) s21_matrix.a -lm $(CHECK_CFLAGS) $(CHECK_LIBS) -fsanitize=address -o test_test
	./test_test
	rm *.a *.out
	rm -rf report

valgrind_test: s21_matrix.a
	gcc $(FLAGS) $(TEST_FILES) s21_matrix.a -lm $(CHECK_CFLAGS) $(CHECK_LIBS) -o test_out
	valgrind --tool=memcheck --leak-check=yes -s ./test_out
	rm *.a test_out
