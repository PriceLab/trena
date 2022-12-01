IMG = bigfimobase
IMG2 = bigfimo
DOCKER = sudo docker

default:
	@echo targets: build.1,2 bash.1,2 build2.noCache testClass
	@echo targets: test.rscript.in.container.small
	@echo targets: test.rscript.in.container.medium
	@echo targets: test.rscript.in.container.huge
build.1:
	$(DOCKER) build -t $(IMG) ./ -f Dockerfile.1

build.2:
	$(DOCKER) build -t $(IMG2) ./ -f Dockerfile.2

build.2.noCache:
	$(DOCKER) build --no-cache -t $(IMG2) ./ -f Dockerfile.2

bash.1:
	$(DOCKER) run -ti --rm $(IMG) bash
bash.2:
	$(DOCKER) run \
           -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results \
           -v /users/pshannon/github/bigFimo/docker/inputs:/usr/local/data/inputs \
           -ti $(IMG2) bash

testClass:
	sudo docker run -ti --rm bigfimo2 R -f /usr/local/scripts/test_BigFimoPackage.R

test.script.in.container.small:
	$(DOCKER) run \
          -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results \
          -v /users/pshannon/github/bigFimo/docker/inputs:/usr/local/data/inputs \
          -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19 4040801 4041801 ZBTB7A 2 1e-5 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results	

# 2 processes, over 1kb, finds 62 hits, ~2kb file
test.rscript.in.container.small:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr19 4040801 4041801 ZBTB7A 2 1e-5 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results								

# 5 processes, 68kb: 715 fimo hits saved in a file ~20k 
test.rscript.in.container.medium:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr1 161165894 161234221 NDUFS2 5 1e-6 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results								

# 20 processes, 1093kb:  2,737,622  fimo hits save in 46,443,187 file, 12 minutes
test.rscript.in.container.huge:
	$(DOCKER) run -v /users/pshannon/github/bigFimo/docker/results:/usr/local/data/results -ti $(IMG2) \
	  Rscript /usr/local/scripts/runBigFimo.R \
              chr1 160653434 161746681 NDUFS2 20 1e-3 \
	      /usr/local/data/human-jaspar2022-hocomoco-core-A.meme \
              /usr/local/data/results								

