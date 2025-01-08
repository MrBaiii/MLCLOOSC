# This Makefile can be used with GNU Make or BSD Make

CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wvla -Wpointer-arith -O3 -march=native -mtune=native
NISTFLAGS += -Wno-unused-result -O3
SOURCES = polyvec.c poly.c ntt.c reduce.c
HEADERS = params.h polyvec.h poly.h ntt.h reduce.h symmetric.h randombytes.h
KECCAK_SOURCES = $(SOURCES) fips202.c symmetric-shake.c
KECCAK_HEADERS = $(HEADERS) fips202.h
AES_SOURCES = $(SOURCES) fips202.c aes256ctr.c symmetric-aes.c
AES_HEADERS = $(HEADERS) fips202.h aes256ctr.h

.PHONY: bench clean

bench_mlcloosc: test/bench.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES)


clean:
	rm -f *~ test/*~ *.gcno *.gcda *.lcov
	rm -f test/bench_mlcloosc
