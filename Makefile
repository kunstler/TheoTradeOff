PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')

all: install

document: roxygen

roxygen: R/RcppExports.R
	@mkdir -p man
	Rscript -e "library(methods); devtools::document()"

R/RcppExports.R: src/CellAuto.cpp
	Rscript -e "Rcpp::compileAttributes(verbose=TRUE)"

install: roxygen
	R CMD INSTALL .

build:
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck

# test:
# 	make -C tests/testthat

.PHONY: attributes document roxygen install clean build check
