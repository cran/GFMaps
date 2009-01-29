\name{GFMap}
\alias{GFMap}
\title{ Visualization technique for interpretation of high-throughput genetic or proteomic experiments. }
\description{
  This technique organizes genes on virtual microarray real-estate according to the prominent biological functionality in the experiment as revealed Gene Set Enrichment analysis.
	   Consistent with enrichment analysis, biological functionality spans Biological Process,Molecular Function and cellular component from Gene Ontology ((GO-BP, GO-MF, GO-CC) as well as KEGG Pathways.
}
\usage{
GFMap(input.ds, class.labels, option, reshuffling.type = "sample.labels", nperm = 100, weighted.score.type = 0, topgs = 10, reverse.sign = FALSE, perm.type = 0, fraction = 1, replace = FALSE, which.GeneRanking = "FC", is.log = FALSE, FUN = mean, bootstrap = TRUE, which.pvalue = "FDR q-val", pvalueCutoff = NULL, use.fast.enrichment.routine = TRUE, importance = "size", relevance = "score", output.directory = "")
}
\arguments{
  \item{input.ds}{ Dataset as an object of class ExpressionSet. }
  \item{class.labels}{ factor level assignment (0 & 1) of datatype numeric }
  \item{option}{ Biological principle for functional organization e.g. GO annotation for Biological Process (BP), Molecular function (MF), Cellular component (CC) or Pathways (KEGG).}
  \item{reshuffling.type}{ Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels") }
  \item{nperm}{ Number of random permutations (default: 100) }
  \item{weighted.score.type}{ Enrichment correlation-based weighting: 0=no weight (KS), 1=standard weigth, 2 = over-weigth (default: 1)  }
  \item{topgs}{ Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10) }
  \item{reverse.sign}{ Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: FALSE) }
  \item{perm.type}{ Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0) }
  \item{fraction}{ Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0) }
  \item{replace}{ Resampling mode (replacement or not replacement). For experts only (default: FALSE) }
  \item{which.GeneRanking}{ 2 metrics for Gene Ranking are available: "FC"- fold change ,"s2n"=Signal 2 noise }
  \item{is.log}{ A logical variable indicating whether the data has been logged.}
  \item{FUN}{ The summary statistics function used to calcuate fold change, the default is set as mean, the user can also use median }
  \item{bootstrap}{ Bootstraping for calculation of normalized enrichmnet score,Default=TRUE }
  \item{which.pvalue}{ 4 p-values: "FDR q-val"-False discovery rate q-values  ,"NOM p-val"-Nominal p-value,"FWER p-val"- Family wise error rate p-values ,"FDR (median)"-FDR q-values from the median of the null distributions.}
  \item{pvalueCutoff}{ A numeric values between zero and one used as a p-value cutoff for p-values (default is NULL) }
  \item{use.fast.enrichment.routine}{ if true it uses a faster version to compute random perm. enrichment "GSEA.EnrichmentScore2" (default is NULL) }
  \item{importance}{ Signifcant categories are categorized on the basis of either 'size' or 'ES' (default is 'size') }
  \item{relevance}{ Assertiveness is communicated on the basis of differential expression "score" (eg FC or s2n) or on the basis of statistical significance p-value 'pvalue' (default is 'size') }
  \item{output.directory}{ Directory where to store output and results (default: .) }
  }
\details{
  This technique organizes genes on virtual microarray real-estate according to the prominent biological functionality in the experiment as revealed by Gene Set Enrichment analysis so as to improve Signal-to-noise ratio and hence, enhance the quality of analysis. 
  The GFMaps is the conceptual analog of a pioneering technique in clinical studies - fMRI . In fMRI, areas of the neuro-anatomy light-up in response to stimulus, revealing the connection between anatomy and cognitive function. Similarly in a GFMaps, parts of the 
  genome light-up to indicate connection with observed phenotypical differential.
}
\value{
  The results of the method are stored in the "output.directory" specified by the user as part of the input parameters. 
  A list containing the 2 dataframes & one image(.JPEG or png)
  \item{Summary }{if bootstrap=TRUE then it contains Gene Set name, Gene set definition, actual gene count for each term tested, Enrichment Score, Normalized Enrichment score (NES) and P value for each category term tested. The rows of the summary table are sorted by the NES values (from maximum positive or negative NES to minimum). Else (bootstrap=FALSE) it will not contain NES, pvalue and rows of the summary table are sorted by Enrichment Score }
  \item{GeneMatrix }{Tt displays the list of genes present in the ExpressionSet along with its details that includes its score (fold change or Signal to noise ratio), P value, Gene Set and Gene Set definition. Each gene is represented by its Entrez Gene ID. 
					Genes that do not map to any Gene Set shows NAs in the corresponding columns}
  \item{Geneomic Function Map(GFMap) }{ An image (.JPEG or pdf) that integrate microarray data with biological annotations in terms of Gene Ontology and pathway information and organize genes in such a way that most assertive genes (i.e high score(FC or S2n) or less pvalue) in biological component are placed at the center of their field, with the less assertive members placed peripherally around this core in successive layers of diminishing assertiveness.}
  }
\references{ 
\itemize{
\item{1.}{ Subramanian, A., Tamayo, P., Mootha, V., K. Mukherjee, S., Ebert, B., L., Gillette, M., A., Paulvich, A., Pomeroy, S., L., Golub, T., R., Lander, E., S., Mesirov, J., P. (2005),Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles \strong{PNAS 102, 15545-15550}}
\item{2.}{ Subramanian, A., Keuhn, H., Gould, J., Tamayo, P., Mesirov, J., P. (2007),GSEA-P: a desktop application for Gene Set Enrichment Analysis \strong{Bioinformatics 23(23), 3251-3253}}
 }
 }
\author{ Dr. Sanjay Bhatikar \email{sanjay.bhatikar@biobase.in}, Kanika Arora \email{kanika.arora@biobase.in} }
\note{ This technique can only handle one factor two level microarray experiment }
\seealso{ R-GSEA \url{http://www.broad.mit.edu/cancer/software/gsea/wiki/index.php/R-GSEA_Readme}}
\examples{

require(ALL)
data(ALL, package = "ALL")
Bcell <- grep("^B", as.character(ALL$BT))
bcrAblOrNegIdx <- which(match(as.character(ALL$mol),c("NEG", "ALL1/AF4"),nomatch=0)!=0)
bcrAblOrNeg <- ALL[, intersect(Bcell, bcrAblOrNegIdx)]
class.labels<-(as.numeric((bcrAblOrNeg$mol.biol!="NEG")))
res<-GFMap(input.ds=bcrAblOrNeg,class.labels=class.labels,nperm=10,option="KEGG",topgs = 20,pvalueCutoff=NULL,which.GeneRanking="s2n",importance="size",output.directory ="")
}

\keyword{ color }
\keyword{ nonparametric }
