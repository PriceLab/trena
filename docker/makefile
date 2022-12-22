IMG = trenabase
IMG2 = trena2
IMG3 = trena
DOCKER = sudo docker
#DOCKER = docker

default:
	@echo targets: build.1,2,3 bash.1,2,3
	@echo targets: test.ft  
	@echo targets: test.trena

build.1:
	$(DOCKER) build -t $(IMG) ./ -f Dockerfile.1

build.2:
	$(DOCKER) build -t $(IMG2) ./ -f Dockerfile.2

build.3:
	$(DOCKER) build -t $(IMG3) ./ -f Dockerfile.3

bash.1:
	$(DOCKER) run -ti --rm $(IMG) bash
bash.3:
	$(DOCKER) run \
           -ti $(IMG3) bash

#            -v /Users/paul/github/trena/docker/results:/usr/local/data/results \
#            -v /Users/paul/github/trena/docker/inputs:/usr/local/data/inputs \

test.ft:
	$(DOCKER) run \
           -v /Users/paul/github/trena/docker/results:/usr/local/data/results \
           -v /Users/paul/github/trena/docker/inputs:/usr/local/data/inputs \
           -ti $(IMG3) R -f /usr/local/scripts/test_FeatureTable.R

test.trena:
	$(DOCKER) run \
           -v /Users/paul/github/trena/docker/results:/usr/local/data/results \
           -v /Users/paul/github/trena/docker/inputs:/usr/local/data/inputs \
           -ti $(IMG3) R -f /usr/local/lib/R/site-library/trena/unitTests/test_EnsembleSolver.R

