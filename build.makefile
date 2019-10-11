quick: roxy install

all:  roxy vig build check

roxy:
	R -e "devtools::document()"
vig:
	R -e "devtools::build_vignettes()"

build:
	(cd ..; R CMD build --no-build-vignettes trena)

install:
	(cd ..; R CMD INSTALL --no-test-load trena)

check:
	(cd ..; R CMD check --no-manual `ls -t trena_* | head -1`)

biocCheck:
	(cd ..; R CMD BiocCheck `ls -t trena_* | head -1`)

unitTests:
	 for x in inst/unitTests/test_*.R; do echo ============== $$x; R -f $$x; done

site:
	R -e "devtools::build_site()"
