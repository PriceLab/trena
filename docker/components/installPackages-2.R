library(devtools)

github.pkgs <- c("Rdatatable/data.table", "PriceLab/ghdb", "PriceLab/trena")
for(pkg in github.pkgs){
    install_github(pkg)
    }
