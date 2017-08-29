#' @title Create a GeneOntologyFilter object
#'
#' @description
#' A GeneOntologyFilter object allows for filtering based on gene ontology (GO) terms. Its
#' associated \code{getCandidates} method uses an organism database and a GO term (e.g. GO:#######)
#' to filter a list of possible regulators factors to those that match the GO term.
#'
#' @include CandidateFilter.R
#' @import methods
#'
#' @name GeneOntologyFilter-class
#' @rdname GeneOntologyFilter-class
#' @aliases GeneOntologyFilter
#----------------------------------------------------------------------------------------------------
.GeneOntologyFilter <- setClass("GeneOntologyFilter",
                            contains="CandidateFilter",
                            slots=c(organismDatabase="OrgDb",
                                    GOTerm="character",
                                    quiet="logical"))
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
# some candidate terms:
#  GO:0010467                        gene expression
#  GO:0097659   nucleic acid-templated transcription
#  GO:0001172           transcription, RNA-templated
#  GO:0006351           transcription, DNA-templated
#----------------------------------------------------------------------------------------------------
GeneOntologyFilter <- function(organismDatabase=org.Hs.eg.db, GOTerm="GO:0006351", quiet=TRUE)
{
    data.file <- system.file(package="trena", "extdata", "human.regulatory.genes.RData")
    load(data.file)
    .GeneOntologyFilter(organismDatabase=organismDatabase, GOTerm=GOTerm, quiet=quiet)

} # GeneOntologyFilter, the constructor
#----------------------------------------------------------------------------------------------------
setMethod("getCandidates", "GeneOntologyFilter",

    function(obj){
       suppressMessages(
          regulatory.genes <- select(org.Hs.eg.db, keys=obj@GOTerm, keytype="GO", columns=c("SYMBOL", "GO"))$SYMBOL
          )
       list(tbl=data.frame(), tfs=regulatory.genes)
       }) # getCandidates

#----------------------------------------------------------------------------------------------------
