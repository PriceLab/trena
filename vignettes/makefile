overview:
	R -e "rmarkdown::render('overview.Rmd')"
	open overview.html

all: overview case_1 case_2 case_3 case_4

case_1:
	R -e "rmarkdown::render('caseStudyOne.Rmd')"
	open caseStudyOne.html

case_2:
	R -e "rmarkdown::render('caseStudyTwo.Rmd')"
	open caseStudyTwo.html

case_3:
	R -e "rmarkdown::render('caseStudyThree.Rmd')"
	open caseStudyThree.html

case_4:
	R -e "rmarkdown::render('caseStudyFour.Rmd')"
	open caseStudyFour.html



vigExtract:
	R -e "library(knitr); purl('overview.Rmd')"

