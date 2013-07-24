include ../../config.mk
include ../../rules.mk

bench_cholesky_ccsd_lambda: $(bindir)/bench-cholesky-ccsd-lambda
$(bindir)/bench-cholesky-ccsd-lambda: cholesky-ccsd-lambda.o

bench_cholesky_ccsd: $(bindir)/bench-cholesky-ccsd
$(bindir)/bench-cholesky-ccsd: cholesky-ccsd.o

bench_cholesky_ccsdt: $(bindir)/bench-cholesky-ccsdt
$(bindir)/bench-cholesky-ccsdt: cholesky-ccsdt.o

bench_ao_ccsd_lambda: $(bindir)/bench-ao-ccsd-lambda
$(bindir)/bench-ao-ccsd-lambda: ao-ccsd-lambda.o
  
bench_ao_ccsd: $(bindir)/bench-ao-ccsd
$(bindir)/bench-ao-ccsd: ao-ccsd.o

bench_ao_ccsdt: $(bindir)/bench-ao-ccsdt
$(bindir)/bench-ao-ccsdt: ao-ccsdt.o
