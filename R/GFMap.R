.packageName <- 'GFMaps'
`GFMap` <-
function(
input.ds, 
class.labels, 
option,
reshuffling.type = "sample.labels", 
nperm = 100, 
weighted.score.type = 0, 
topgs = 10, 
reverse.sign = FALSE,  
perm.type = 0, 
fraction = 1.0, 
replace = FALSE,
which.GeneRanking="FC",
is.log=FALSE,
FUN=mean,
bootstrap=TRUE,
which.pvalue="FDR q-val",
pvalueCutoff=NULL,
use.fast.enrichment.routine = TRUE ,
importance="size",
relevance="score",
output.directory ="") {

# This is a methodology for the analysis of global molecular profiles called Gene Set Enrichment Analysis (GSEA). It determines 
# whether an a priori defined set of genes shows statistically significant, concordant differences between two biological 
# states (e.g. phenotypes). GSEA operates on all genes from an experiment, rank ordered by the signal to noise ratio and 
# determines whether members of an a priori defined gene set are nonrandomly distributed towards the top or bottom of the 
# list and thus may correspond to an important biological process. To assess significance the program uses an empirical 
# permutation procedure to test deviation from random that preserves correlations between genes. 
#
# 
#
# Inputs:
#   input.ds: Expressionset 
#   class.labels: factor level assignment (0 & 1) of datatype numeric
#   option=GO ("BP","CC","MF")
#   reshuffling.type: Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels") 
#   nperm: Number of random permutations (default: 100) 
#   weighted.score.type: Enrichment correlation-based weighting: 0=no weight (KS), 1=standard weigth, 2 = over-weigth (default: 1) 
#   topgs: Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10) 
#   reverse.sign: Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: FALSE) 
#   perm.type: Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0) 
#   fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0) 
#   replace: Resampling mode (replacement or not replacement). For experts only (default: FALSE)
#   GeneRanking: 2 metrics for Gene Ranking are available: "FC"- fold change ,"s2n"=Signal 2 noise
#   is.log:A logical variable indicating whether the data has been logged.
#   FUN:The summary statistics function used to calcuate fold change, the default is set as mean, the user can also use median
#   bootstrap:Bootstraping for calculation of normalized enrichmnet score,Default=TRUE
#   which.pvalue: 4 p-values: "FDR q-val"-False discovery rate q-values  ,"NOM p-val"-Nominal p-value,"FWER p-val"- Family wise error rate p-values ,"FDR (median)"-FDR q-values from the median of the null distributions.
#   pvalueCutoff:A numeric values between zero and one used as a p-value cutoff for p-values (default is NULL)
#   use.fast.enrichment.routine: if true it uses a faster version to compute random perm. enrichment "GSEA.EnrichmentScore2"  
#   importance: Signifcant categories are categorized on the basis of either 'size' or 'ES'
#   relevance: Assertiveness is communicated on the basis of differential expression "score" (eg FC or s2n) or on the basis of statistical significance p-value 'pvalue' (default is 'size')
#   output.directory: Directory where to store output and results (default: .) 

#Output:
#Output:
#  The results of the method are stored in the "output.directory" specified by the user as part of the input parameters. 
#  The results files are:
#A list containing the 2 dataframes & one image(png)
# a) Summary :It conatins following components:
# The format (columns) for the global result files is as follows.
# GS : Gene set name.
# GS.term: Gene set term name
# SIZE : Size of the set in genes.
# ES : Enrichment score.
# NES : Normalized (multiplicative rescaling) normalized enrichment score.
# p-val : p-values for each category term tested.
# The rows are sorted by the NES values (from maximum positive or negative NES to minimum)


# b) GeneMatrix:It conatins following components:
# Probeset_ids:All the Probesetis in the Expressionset(input)
# Term_ID: gene set Id associated with Id.
# GeneScore:either signal 2 noise or Fold change (depending on user specification) for each probesetids
# Genepvalue:pvalue associated with score(calculated via nonparametrized bootstraping)
# Term:Name of the gene set.

 

#----------------------------------------------------------------------------------------------------
#  loading libraries
#----------------------------------------------------------------------------------------------------
  
  library(affy)
  library(topGO)
  affylib<-annotation(input.ds)
  lib_db<-paste(affylib,".db",sep="")
  library(package=lib_db,character.only=TRUE)
  library(KEGG.db)
  library(plotrix)

  print(" *** Running GSEA Analysis...")
   
   

#---------------------------------------------------------------------------------------------  
# Start of GSEA methodology 
#---------------------------------------------------------------------------------------------
  dataset<-data.frame(exprs(input.ds))
  gene.labels <- row.names(dataset)
  sample.names <- names(dataset)
  A <- data.matrix(dataset)
  cols <- length(A[1,])
  rows <- length(A[,1])

  # sort samples according to phenotype
 
  col.index <- order(class.labels, decreasing=FALSE)
  class.labels <- class.labels[col.index]
  sample.names <- sample.names[col.index]
  for (j in 1:rows) {
     A[j, ] <- A[j, col.index]
  }
 colnames(A) <- sample.names
 N <- length(A[,1])
 Ns <- length(A[1,])

 
 #------------------------------------------------------------------------------
 # 1)  Function definition of function Geneset to do GO or KEGG analysis
 #-------------------------------------------------------------------------------
 Geneset<-function(option,affylib)
 {
    if(option=="KEGG")
      {
        path<-paste(affylib,"PATH2PROBE",sep="")
	    gs<-eval(parse(text=paste("gs<-as.list(",path,")")))
	    gs.terms<-as.character((mget(names(gs),envir=KEGGPATHID2NAME)))
	} else {
        go2<-paste(affylib,"GO2PROBE",sep="")
	    GO2gene<-eval(parse(text=paste("GO2gene<-as.list(",go2,")")))
        gs<- annFUN.GO2genes(whichOnto = option, GO2gene = GO2gene)
	    gs.terms<-as.character(lapply(mget(names(gs),envir=GOTERM),Term))
	 }
	 return(list(gs=gs,gs.terms=gs.terms))
}
 
#-------------------------------------------------------------------------------------------------------------------- 
# Read input gene set database(Gene set list for enrichment function
#-------------------------------------------------------------------------------------------------------------------
  G<-Geneset(option,affylib)    # call to function Geneset
  gs<-G$gs
  gs.terms<-G$gs.terms
  Ng<-length(gs)
  gs.names<-names(gs)
  size.G <-as.numeric(lapply(gs,length))
 
#----------------------------------------------------------------------------------------------------------------------
# 2)   GeneRanking (ranked List for enrichemnt function)
#---------------------------------------------------------------------------------------------------------------------- 
# This function calls two function :: 
# 2a)  GeneRanking.FC 
# 2b)  GSEA.GeneRanking.s2n


#---------------------------------------------------------------------------------------
# 2a)	GSEA.GeneRanking.FC ( Fold change)
#---------------------------------------------------------------------------------------
# This function calls 1 function:: Compute.FC
# 2aa) Function Defintion of Compute.FC.R
#---------------------------------------------------------------------------------------


Compute.FC<-function(X,L = NULL, is.log = FALSE, FUN = mean)
{
         if (is.vector(X))
            X <- matrix(X, byrow = TRUE)
        if (is.null(L))
            L <- rep(0, ncol(X))
        nL <-length(unique(L))
        nr <-nrow(X)
        nc <- ncol(X)
        if (nL == 1) {
            fc <- apply(X, 1, FUN, na.rm = TRUE)
            return(fc)
        }
        else if (nL == 2) {
            G1 <- X[, L == 0]
            G2 <- X[, L == 1]
            m1 <- apply(G1, 1, FUN, na.rm = TRUE)
            m2 <- apply(G2, 1, FUN, na.rm = TRUE)
            if (is.log)
                fc <- m1 - m2
            else fc <- log2(m1/m2)
            fc
        }
        else {
            m <- t(apply(X, 1, function(z) {
                tapply(z, L, FUN)
            }))
            if (is.log)
                fc <- apply(m, 1, function(z) {
                  max(z) - min(z)
                })
            else fc <- apply(m, 1, function(z) {
                log2(max(z)/min(z))
            })
            return(fc)
        }
    }

#------------------------------------------------------------------------------------------
# 2a)	Function definition of function GSEA.GeneRanking.FC ( Fold change)
#---------------------------------------------------------------------------------------

GSEA.GeneRanking.FC <- function(A, class.labels,nperm, is.log=FALSE,FUN=mean,reverse.sign=FALSE) { 

# This function ranks the genes according to the fold change.
#
# Inputs:
#   A: Matrix of gene expression values (rows are genes, columns are samples) 
#   class.labels: Phenotype of class disticntion of interest. A vector of binary labels having first the 1's and then the 0's 
#   nperm: Number of random permutations/bootstraps to perform 
#   is.log=A logical variable indicating whether the data has been logged.
#   FUN=The summary statistics function used to calcuate fold change, the default is set as mean, the user can also use median
#
# Outputs:
#   FC.matrix: Matrix with random permuted or bootstraps Fold change values(rows are genes, columns are permutations or bootstrap subsamplings
#   obs.FC.matrix: Matrix with observed Fold change
#   order.matrix: Matrix with the orderings that will sort the columns of the obs.FC.matrix in decreasing FC order
#   obs.order.matrix: Matrix with the orderings that will sort the columns of the FC.matrix in decreasing FC order
#


     A <- A + 0.00000001

     N <- length(A[,1])
     Ns <- length(A[1,])

     


     order.matrix <- matrix(0, nrow = N, ncol = nperm)
     obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
     FC.matrix <- matrix(0, nrow = N, ncol = nperm)
     obs.FC.matrix <- matrix(0, nrow = N, ncol = nperm)

     
     fc_val<-Compute.FC(A,L=class.labels,is.log=is.log,FUN=FUN)
     for (r in 1:nperm) {
        reshuffled.class.labels<- sample(class.labels)
		# compute FC for the random permutation matrix
        FC.matrix[,r]<-Compute.FC(A,L=reshuffled.class.labels,is.log=is.log,FUN=FUN)
		# compute FC for the "observed" permutation matrix
		obs.FC.matrix[,r]<-fc_val
		}

     
     if (reverse.sign == TRUE) {
        FC.matrix <- - FC.matrix
     }
     gc()

     for (r in 1:nperm) {
        order.matrix[, r] <- order(FC.matrix[, r], decreasing=TRUE)            
     }


     for (r in 1:nperm) {
        obs.order.matrix[,r] <- order(obs.FC.matrix[,r], decreasing=TRUE)            
     }
    
	

     return(list(FC.matrix = FC.matrix, 
                 obs.FC.matrix = obs.FC.matrix, 
                 order.matrix = order.matrix,
                 obs.order.matrix = obs.order.matrix))
}
#-------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------
# 2b)	Function definition of function GeneRanking.s2n ( Signal 2 Noise Ratio)
#-------------------------------------------------------------------------------------------
GSEA.GeneRanking.s2n<- function(A, class.labels, nperm, permutation.type = 0, fraction=1.0, replace=FALSE, reverse.sign= FALSE) { 

# This function ranks the genes according to the signal to noise ratio for the actual phenotype and also random permutations and bootstrap  
# subsamples of both the observed and random phenotypes. It uses matrix operations to implement the signal to noise calculation 
# in stages and achieves fast execution speed. It supports two types of permutations: random (unbalanced) and balanced. 
# It also supports subsampling and bootstrap by using masking and multiple-count variables.  When "fraction" is set to 1 (default)
# the there is no subsampling or boostrapping and the matrix of observed signal to noise ratios will have the same value for 
# all permutations. This is wasteful but allows to support all the multiple options with the same code. Notice that the second 
# matrix for the null distribution will still have the values for the random permutations 
# (null distribution). This mode (fraction = 1.0) is the defaults, the recommended one and the one used in the examples.
# It is also the one that has be tested more thoroughly. The resampling and boostrapping options are intersting to obtain 
# smooth estimates of the observed distribution but its is left for the expert user who may want to perform some sanity 
# checks before trusting the code.
#
# Inputs:
#   A: Matrix of gene expression values (rows are genes, columns are samples) 
#   class.labels: Phenotype of class disticntion of interest. A vector of binary labels having first the 1's and then the 0's 
#   nperm: Number of random permutations/bootstraps to perform 
#   permutation.type: Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0) 
#   fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0) 
#   replace: Resampling mode (replacement or not replacement). For experts only (default: FALSE) 
#   reverse.sign: Reverse direction of gene list (default = FALSE)
#
# Outputs:
#   s2n.matrix: Matrix with random permuted or bootstraps signal to noise ratios (rows are genes, columns are permutations or bootstrap subsamplings
#   obs.s2n.matrix: Matrix with observed signal to noise ratios (rows are genes, columns are boostraps subsamplings. If fraction is set to 1.0 then all the columns have the same values
#   order.matrix: Matrix with the orderings that will sort the columns of the obs.s2n.matrix in decreasing s2n order
#   obs.order.matrix: Matrix with the orderings that will sort the columns of the s2n.matrix in decreasing s2n order
#

     A <- A + 0.00000001

     N <- length(A[,1])
     Ns <- length(A[1,])

     subset.mask <- matrix(0, nrow=Ns, ncol=nperm)
     reshuffled.class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
     reshuffled.class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
     class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
     class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
	
     order.matrix <- matrix(0, nrow = N, ncol = nperm)
     obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
     s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
     obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)


     M1 <- matrix(0, nrow = N, ncol = nperm)
     M2 <- matrix(0, nrow = N, ncol = nperm)
     S1 <- matrix(0, nrow = N, ncol = nperm)
     S2 <- matrix(0, nrow = N, ncol = nperm)

	 
     gc()


     class1.size <- sum(class.labels==0)
     class2.size <- sum(class.labels==1)
     class1.index <- seq(1, class1.size, 1)
     class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)

     for (r in 1:nperm) {
        class1.subset <- sample(class1.index, size = ceiling(class1.size*fraction), replace = replace)
        class2.subset <- sample(class2.index, size = ceiling(class2.size*fraction), replace = replace)
        class1.subset.size <- length(class1.subset)
        class2.subset.size <- length(class2.subset)
        subset.class1 <- rep(0, class1.size)
        subset.class1[(class1.index %in% class1.subset)]<-1
        subset.class2 <- rep(0, class2.size)
        subset.class2[(class2.index %in% class2.subset)]<-1
        subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
        fraction.class1 <- class1.size/Ns
        fraction.class2 <- class2.size/Ns

        if (permutation.type == 0) { # random (unbalanced) permutation
           full.subset <- c(class1.subset, class2.subset)
           label1.subset <- sample(full.subset, size = Ns * fraction.class1)
           reshuffled.class.labels1[(1:Ns %in% label1.subset),r]<-1
		   reshuffled.class.labels2[(1:Ns%in% full.subset[which(!(full.subset %in% label1.subset))]),r]<-1
		   class.labels1[((1:Ns %in% full.subset)[1:class1.size]),r]<-1
           class.labels1[(class1.size+1):Ns,r]<-0
           class.labels2[(1:Ns %in% full.subset),r]<-1
		   class.labels2[1:class1.size,r]<-0
			   
        } else if (permutation.type == 1) { # proportional (balanced) permutation

           class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size*fraction.class1))
           class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size*fraction.class1))
           for (i in 1:Ns) {
               if (i <= class1.size) {
				  m1<-sum(class1.label1.subset %in% i)
                  m2 <- sum(class1.subset %in% i)
                  reshuffled.class.labels1[i, r] <- m1
                  reshuffled.class.labels2[i, r] <- m2 - m1
                  class.labels1[i, r] <- m2
                  class.labels2[i, r] <- 0
               } else {
                  m1 <- sum(class2.label1.subset %in% i)
                  m2 <- sum(class2.subset %in% i)
                  reshuffled.class.labels1[i, r] <- m1
                  reshuffled.class.labels2[i, r] <- m2 - m1
                  class.labels1[i, r] <- 0
                  class.labels2[i, r] <- m2
               }
           }
        }
    }

# compute S2N for the random permutation matrix
     
     P <- reshuffled.class.labels1 * subset.mask
     n1 <- sum(P[,1])         
     M1 <- A %*% P
     M1 <- M1/n1      
     gc()
     A2 <- A*A        
     S1 <- A2 %*% P   
     S1 <- S1/n1 - M1*M1    
     S1 <- sqrt(abs((n1/(n1-1)) * S1))   
     gc()
     P <- reshuffled.class.labels2 * subset.mask
     n2 <- sum(P[,1])           
     M2 <- A %*% P           
     M2 <- M2/n2          
     gc()
     A2 <- A*A           
     S2 <- A2 %*% P      
     S2 <- S2/n2 - M2*M2 
     S2 <- sqrt(abs((n2/(n2-1)) * S2))
     rm(P)
     rm(A2)
     gc()

       # small sigma "fix" as used in GeneCluster
         S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
         S2 <- ifelse(S2 == 0, 0.2, S2)
         S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
         S1 <- ifelse(S1 == 0, 0.2, S1)
         gc()
     

     M1 <- M1 - M2
     rm(M2)
     gc()
     S1 <- S1 + S2
     rm(S2)
     gc()

     s2n.matrix <- M1/S1

     if (reverse.sign == TRUE) {
        s2n.matrix <- - s2n.matrix
     }
     gc()

     for (r in 1:nperm) {
        order.matrix[, r] <- order(s2n.matrix[, r], decreasing=TRUE)            
     }

# compute S2N for the "observed" permutation matrix

     P <- class.labels1 * subset.mask
     n1 <- sum(P[,1])         
     M1 <- A %*% P
     M1 <- M1/n1      
     gc()
     A2 <- A*A        
     S1 <- A2 %*% P   
     S1 <- S1/n1 - M1*M1    
     S1 <- sqrt(abs((n1/(n1-1)) * S1))   
     gc()
     P <- class.labels2 * subset.mask
     n2 <- sum(P[,1])           
     M2 <- A %*% P           
     M2 <- M2/n2          
     gc()
     A2 <- A*A           
     S2 <- A2 %*% P      
     S2 <- S2/n2 - M2*M2 
     S2 <- sqrt(abs((n2/(n2-1)) * S2))
     rm(P)
     rm(A2)
     gc()

      # small sigma "fix" as used in GeneCluster
         S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
         S2 <- ifelse(S2 == 0, 0.2, S2)
         S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
         S1 <- ifelse(S1 == 0, 0.2, S1)
         gc()

     M1 <- M1 - M2
     rm(M2)
     gc()
     S1 <- S1 + S2
     rm(S2)
     gc()

     obs.s2n.matrix <- M1/S1
     gc()

     if (reverse.sign == TRUE) {
        obs.s2n.matrix <- - obs.s2n.matrix
     }

     for (r in 1:nperm) {
        obs.order.matrix[,r] <- order(obs.s2n.matrix[,r], decreasing=TRUE)            
     }

     return(list(s2n.matrix = s2n.matrix, 
                 obs.s2n.matrix = obs.s2n.matrix, 
                 order.matrix = order.matrix,
                 obs.order.matrix = obs.order.matrix))
}


#----------------------------------------------------------------------------------------------------------------------
# 2)  call GeneRanking function (ranked List for enrichemnt function)
#---------------------------------------------------------------------------------------------------------------------- 
GeneRanking<-function(A,class.labels,which.GeneRanking,nperm,is.log,FUN,reverse.sign,replace=replace,perm.type, fraction=fraction)
	   {
	        N<-(dim(A)[1])
		    correl.matrix <- matrix(nrow = N, ncol = nperm)
            obs.correl.matrix <- matrix(nrow = N, ncol = nperm)
            order.matrix <- matrix(nrow = N, ncol = nperm)
            obs.order.matrix <- matrix(nrow = N, ncol = nperm)
			obs.val <- vector(length=N, mode="numeric")
	        #pval.gene<- vector(length=N, mode="numeric")
		    if(which.GeneRanking=="FC")
		      {
				nperm.per.call <- 100
                n.groups <- nperm %/% nperm.per.call
                n.rem <- nperm %% nperm.per.call
                n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
                n.ends <- cumsum(n.perms)
                n.starts <- n.ends - n.perms + 1

                if (n.rem == 0) {
                  n.tot <- n.groups
                } else {
                  n.tot <- n.groups + 1
	              }

             for (nk in 1:n.tot) {
                 call.nperm <- n.perms[nk]

                 print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", n.starts[nk], "--", n.ends[nk], sep=" "))
                 O<-GSEA.GeneRanking.FC(A, class.labels,call.nperm,is.log=is.log,FUN=FUN,reverse.sign = reverse.sign)
		         gc()
		         order.matrix[,n.starts[nk]:n.ends[nk]] <- O$order.matrix
                 obs.order.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
                 correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$FC.matrix
                 obs.correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.FC.matrix
                 rm(O)
                }					 
				 obs.val<- apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
				 obs.index <- order(obs.val, decreasing=TRUE)            
				 obs.val  <- sort(obs.val, decreasing=TRUE)            

              for (r in 1:nperm) {
				 correl.matrix[, r] <- correl.matrix[order.matrix[,r], r]
                 }
              for (r in 1:nperm) {
                 obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[,r], r]
                 }
		   
		      return(list(obs.val=obs.val,obs.index=obs.index,order.matrix=order.matrix,obs.order.matrix=obs.order.matrix,correl.matrix=correl.matrix,
		      obs.correl.matrix=obs.correl.matrix))
		   
		  } else {
		  
	             nperm.per.call <- 100
                 n.groups <- nperm %/% nperm.per.call
                 n.rem <- nperm %% nperm.per.call
                 n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
                 n.ends <- cumsum(n.perms)
                 n.starts <- n.ends - n.perms + 1

                 if (n.rem == 0) {
                   n.tot <- n.groups
                 } else {
                   n.tot <- n.groups + 1
	              }
	
		      for (nk in 1:n.tot) {
                 call.nperm <- n.perms[nk]

                 print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", n.starts[nk], "--", n.ends[nk], sep=" "))
  
                 O <- GSEA.GeneRanking.s2n(A, class.labels,call.nperm, permutation.type = perm.type,fraction=fraction, replace=replace, reverse.sign = reverse.sign)
                 gc()

                 order.matrix[,n.starts[nk]:n.ends[nk]] <- O$order.matrix
                 obs.order.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
                 correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$s2n.matrix
                 obs.correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.s2n.matrix
                   rm(O)
                  }
		         obs.s2n <- apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
                 obs.index <- order(obs.s2n, decreasing=TRUE)            
                 obs.s2n   <- sort(obs.s2n, decreasing=TRUE)            

              for (r in 1:nperm) {
                   correl.matrix[, r] <- correl.matrix[order.matrix[,r], r]
                 }
             for (r in 1:nperm) {
                  obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[,r], r]
                 }
		        
		     return(list(obs.val=obs.s2n,obs.index=obs.index,order.matrix=order.matrix,obs.order.matrix=obs.order.matrix,correl.matrix=correl.matrix,
		     obs.correl.matrix=obs.correl.matrix))
		   
		}
		
		
}

#------------------------------------------------------------------------------------
# Function defintion of function pvalue
#------------------------------------------------------------------------------------
pvalue<-function(correl.matrix,obs.val,nperm)
{
  pval.gene<- vector(length=length(obs.val), mode="numeric")
  for(i in 1:length(obs.val))
    {
	  pval.gene[i]<-(sum(correl.matrix[i,] > obs.val[i]))/nperm
	}
	return(pval.gene=pval.gene)
 }

		
  obs.val <- vector(length=N, mode="numeric")
  correl.matrix <- matrix(nrow = N, ncol = nperm)
  obs.correl.matrix <- matrix(nrow = N, ncol = nperm)
  order.matrix <- matrix(nrow = N, ncol = nperm)
  obs.order.matrix <- matrix(nrow = N, ncol = nperm)
  pval.gene<- vector(length=N, mode="numeric")
	 
  O<-GeneRanking(A,class.labels,which.GeneRanking,nperm,is.log,FUN,reverse.sign,replace=replace,perm.type,fraction=fraction)
  correl.matrix<-O$correl.matrix
  obs.correl.matrix<-O$obs.correl.matrix
  order.matrix<-O$order.matrix
  obs.order.matrix<-O$obs.order.matrix
  obs.val<-O$obs.val
  obs.index<-O$obs.index 
  pval.gene<-pvalue(correl.matrix,obs.val,nperm)
  gc()
  rm(O);
	 
	
#------------------------------------------------------------------------------------------------------------------
# Enrichemnt score calculation
# This calls function GSEA.EnrichmentScore
#-------------------------------------------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------------------
# 3)	Function definition of function GSEA.EnrichmentScore 
# 3a)   Function definition of GSEA.EnrichmentScore
# 3b)	Function definition of GSEA.EnrichmentScore2
#-------------------------------------------------------------------------------------------
GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
#
# Computes the weighted GSEA score of gene.set in gene.list. 
# The weighted score type is the exponent of the correlation 
# weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
# necessary to input the correlation vector with the values in the same order as in the gene list.
#
# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
#   correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
#   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
#   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
#


   tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
   no.tag.indicator <- 1 - tag.indicator 
   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 
   if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
   }
   alpha <- weighted.score.type
   correl.vector <- abs(correl.vector**alpha)
   sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
   norm.tag    <- 1.0/sum.correl.tag
   norm.no.tag <- 1.0/Nm
   RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
   max.ES <- max(RES)
   min.ES <- min(RES)
   if (max.ES > - min.ES) {
#      ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
   } else {
#      ES <- min.ES
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
   }
   return(list(ES = ES, arg.ES = arg.ES))    
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
#
# Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
# GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
# This call is intended to be used to asses the enrichment of random permutations rather than the 
# observed one.
# The weighted score type is the exponent of the correlation 
# weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
# necessary to input the correlation vector with the values in the same order as in the gene list.
#
# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
#  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#


   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 

   loc.vector <- vector(length=N, mode="numeric")
   peak.res.vector <- vector(length=Nh, mode="numeric")
   valley.res.vector <- vector(length=Nh, mode="numeric")
   tag.correl.vector <- vector(length=Nh, mode="numeric")
   tag.diff.vector <- vector(length=Nh, mode="numeric")
   tag.loc.vector <- vector(length=Nh, mode="numeric")

   loc.vector[gene.list] <- seq(1, N)
   tag.loc.vector <- loc.vector[gene.set]

   tag.loc.vector <- sort(tag.loc.vector, decreasing = FALSE)

   if (weighted.score.type == 0) {
      tag.correl.vector <- rep(1, Nh)
   } else if (weighted.score.type == 1) {
      tag.correl.vector <- correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
   } else if (weighted.score.type == 2) {
      tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
   } else {
      tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
      tag.correl.vector <- abs(tag.correl.vector)
   }

   norm.tag <- 1.0/sum(tag.correl.vector)
   tag.correl.vector <- tag.correl.vector * norm.tag
   norm.no.tag <- 1.0/Nm
   if(Nh==1)
   {
   tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
   } else {
     tag.diff.vector[1] <- (tag.loc.vector[1] - 1)  
     tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
   }
   tag.diff.vector <- tag.diff.vector * norm.no.tag
   peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
   valley.res.vector <- peak.res.vector - tag.correl.vector
   max.ES <- max(peak.res.vector)
   min.ES <- min(valley.res.vector)
   ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)

   return(list(ES = ES))

}

#--------------------------------------------------------------------------------------------

  #Obs.indicator <- matrix(nrow= Ng, ncol=N)
  #Obs.RES <- matrix(nrow= Ng, ncol=N)
  Obs.ES <- vector(length = Ng, mode = "numeric")
  Obs.arg.ES <- vector(length = Ng, mode = "numeric")	  
   
  gene.list2 <- obs.index  
	  
  for (i in 1:Ng) {
    print(paste("Computing observed enrichment for gene set:", i, gs.names[i], sep=" ")) 
    gene.set <- gs[[i]]
    gene.set2 <- vector(length=length(gene.set), mode = "numeric")
    gene.set2 <- match(gene.set, gene.labels)

    GSEA.results <- GSEA.EnrichmentScore(gene.list=gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector = obs.val)
       
    Obs.ES[i] <- GSEA.results$ES
    Obs.arg.ES[i] <- GSEA.results$arg.ES
    #Obs.RES[i,] <- GSEA.results$RES
    #Obs.indicator[i,] <- GSEA.results$indicator
       
     }
#------------------------------------------------------------------------------------------------------------------
# 4)	BOOTSTRAP
# Bootstrap_NES calls::
# a)	Which_p.value
#-------------------------------------------------------------------------------------------------------------------	 
 	 
#	 Which_p.value
#-------------------------------------------------------------------------------------------------
	   
  Nominal_pvalue<-function(Ng,phi,nperm,Obs.ES)
    {
        
	
	     # Find nominal p-values       

		 print("Computing nominal p-values...")

         p.vals <- matrix(0, nrow = Ng, ncol = 1)
	 
	     for (i in 1:Ng) {
          pos.phi <- NULL
          neg.phi <- NULL
            for (j in 1:nperm) {
              if (phi[i, j] >= 0) {
               pos.phi <- c(pos.phi, phi[i, j]) 
             } else {
              neg.phi <- c(neg.phi, phi[i, j]) 
             }
           }
         ES.value <- Obs.ES[i]
         if (ES.value >= 0) {
          p.vals[i] <- signif(sum(pos.phi >= ES.value)/length(pos.phi), digits=5)
          } else {
         p.vals[i] <- signif(sum(neg.phi <= ES.value)/length(neg.phi), digits=5)
         }
       }
	   return(p.vals)
	}
	
	
	FWER<-function(Ng,nperm,phi.norm,Obs.ES.norm)
	 {
	
	     # Compute FWER p-vals

         print("Computing FWER p-values...")

         p.vals <- matrix(0, nrow = Ng, ncol = 1)
         max.ES.vals.p <- NULL
         max.ES.vals.n <- NULL
         for (j in 1:nperm) {
         pos.phi <- NULL
         neg.phi <- NULL
         for (i in 1:Ng) {
          if (phi.norm[i, j] >= 0) {
             pos.phi <- c(pos.phi, phi.norm[i, j]) 
          } else {
             neg.phi <- c(neg.phi, phi.norm[i, j]) 
          }
         }
          if (length(pos.phi) > 0) {
          max.ES.vals.p <- c(max.ES.vals.p, max(pos.phi))
         }
         if (length(neg.phi) > 0) {
          max.ES.vals.n <- c(max.ES.vals.n, min(neg.phi))
         }
       }
       for (i in 1:Ng) {
       ES.value <- Obs.ES.norm[i]
       if (Obs.ES.norm[i] >= 0) {
          p.vals[i] <- signif(sum(max.ES.vals.p >= ES.value)/length(max.ES.vals.p), digits=5)
         } else {
          p.vals[i] <- signif(sum(max.ES.vals.n <= ES.value)/length(max.ES.vals.n), digits=5)
         }
       }
	   return(p.vals)
	 }

	 
	 FDR.mean<-function(Ng,Obs.ES.norm,gs.names,phi.norm,obs.phi.norm,nperm)
	 {
	      # Compute FDRs 

          print("Computing FDR q-values...")

         NES <- vector(length=Ng, mode="numeric")
         phi.norm.mean  <- vector(length=Ng, mode="numeric")
         obs.phi.norm.mean  <- vector(length=Ng, mode="numeric")
         phi.norm.mean  <- vector(length=Ng, mode="numeric")
         obs.phi.mean  <- vector(length=Ng, mode="numeric")
         FDR.mean <- vector(length=Ng, mode="numeric")
      
     

         Obs.ES.index <- order(Obs.ES.norm, decreasing=TRUE)
         Orig.index <- seq(1, Ng)
         Orig.index <- Orig.index[Obs.ES.index]
         Orig.index <- order(Orig.index, decreasing=FALSE)
         Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
         gs.names.sorted <- gs.names[Obs.ES.index]

          for (k in 1:Ng) {
            NES[k] <- Obs.ES.norm.sorted[k]
            ES.value <- NES[k]
            count.col <- vector(length=nperm, mode="numeric")
            obs.count.col <- vector(length=nperm, mode="numeric")
            for (i in 1:nperm) {
              phi.vec <- phi.norm[,i]
              obs.phi.vec <- obs.phi.norm[,i]
            if (ES.value >= 0) {
               count.col.norm <- sum(phi.vec >= 0)
               obs.count.col.norm <- sum(obs.phi.vec >= 0)
               count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 0)
               obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0)
            } else {
               count.col.norm <- sum(phi.vec < 0)
               obs.count.col.norm <- sum(obs.phi.vec < 0)
               count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 0)
               obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
            }
         }
        phi.norm.mean[k] <- mean(count.col)
        obs.phi.norm.mean[k] <- mean(obs.count.col)
        FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
      }

        #adjust q-values

         pos.nes <- length(NES[NES >= 0])
         min.FDR.mean <- FDR.mean[pos.nes]
         for (k in seq(pos.nes - 1, 1, -1)) {
              if (FDR.mean[k] < min.FDR.mean) {
                  min.FDR.mean <- FDR.mean[k]
              }
              if (min.FDR.mean < FDR.mean[k]) {
                  FDR.mean[k] <- min.FDR.mean
              }
         }

         neg.nes <- pos.nes + 1
         min.FDR.mean <- FDR.mean[neg.nes]
         for (k in seq(neg.nes + 1, Ng)) {
             if (FDR.mean[k] < min.FDR.mean) {
                 min.FDR.mean <- FDR.mean[k]
             }
             if (min.FDR.mean < FDR.mean[k]) {
                 FDR.mean[k] <- min.FDR.mean
             }
         }
    

      obs.phi.norm.mean.sorted <- obs.phi.norm.mean[Orig.index]
      phi.norm.mean.sorted <- phi.norm.mean[Orig.index]
      FDR.mean.sorted <- FDR.mean[Orig.index]
      return(FDR.mean.sorted)
	
}



   FDR.median<-function(Ng,Obs.ES.norm,Obs.ES.index,gs.names,phi.norm,obs.phi.norm,nperm)
   {
	  # Compute FDRs 

      print("Computing FDR q-values...")

      NES <- vector(length=Ng, mode="numeric")
      phi.norm.median  <- vector(length=Ng, mode="numeric")
      obs.phi.norm.median  <- vector(length=Ng, mode="numeric")
      FDR.median <- vector(length=Ng, mode="numeric")
      phi.norm.median.d <- vector(length=Ng, mode="numeric")
      obs.phi.norm.median.d <- vector(length=Ng, mode="numeric")

      Obs.ES.index <- order(Obs.ES.norm, decreasing=TRUE)
      Orig.index <- seq(1, Ng)
      Orig.index <- Orig.index[Obs.ES.index]
      Orig.index <- order(Orig.index, decreasing=FALSE)
      Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
      gs.names.sorted <- gs.names[Obs.ES.index]

      for (k in 1:Ng) {
         NES[k] <- Obs.ES.norm.sorted[k]
         ES.value <- NES[k]
         count.col <- vector(length=nperm, mode="numeric")
         obs.count.col <- vector(length=nperm, mode="numeric")
         for (i in 1:nperm) {
            phi.vec <- phi.norm[,i]
            obs.phi.vec <- obs.phi.norm[,i]
            if (ES.value >= 0) {
               count.col.norm <- sum(phi.vec >= 0)
               obs.count.col.norm <- sum(obs.phi.vec >= 0)
               count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 0)
               obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0)
            } else {
               count.col.norm <- sum(phi.vec < 0)
               obs.count.col.norm <- sum(obs.phi.vec < 0)
               count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 0)
               obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
            }
         }
        phi.norm.median[k] <- median(count.col)
        obs.phi.norm.median[k] <- median(obs.count.col)
        FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 1, phi.norm.median[k]/obs.phi.norm.median[k], 1)
      }

    # adjust q-values

      
         pos.nes <- length(NES[NES >= 0])
         min.FDR.median <- FDR.median[pos.nes]
         neg.nes <- pos.nes + 1
         min.FDR.median <- FDR.median[neg.nes]
         FDR.median.sorted <- FDR.median[Orig.index]
		 return(FDR.median.sorted)
    }
	
#-------------------------------------------------------------------------------------------------
# Bootstrap_NES
#-------------------------------------------------------------------------------------------------

Bootstrap_NES<-function(nperm,Ng,reshuffling.type,gs,gene.labels,gs.names,order.matrix,obs.order.matrix,correl.matrix,obs.correl.matrix,gene.list2,use.fast.enrichment.routine,weighted.score.type,fraction,Obs.ES,which.pvalue)
{
      # Compute enrichment for random permutations 
      phi <- matrix(nrow = Ng, ncol = nperm)
      phi.norm <- matrix(nrow = Ng, ncol = nperm)
      obs.phi <- matrix(nrow = Ng, ncol = nperm)
	  obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)
	  Obs.ES.norm <- vector(length = Ng, mode = "numeric")
	  
	  
      if (reshuffling.type == "sample.labels") { # reshuffling phenotype labels
      for (i in 1:Ng) {
        print(paste("Computing random permutations' enrichment for gene set:", i, gs.names[i], sep=" ")) 
        gene.set <- gs[[i]]
        #gene.set2 <- vector(length=length(gene.set), mode = "numeric")
		gene.set2<-NULL
        gene.set2 <- which(gene.labels %in% gene.set)
        for (r in 1:nperm) {
            gene.list2 <- order.matrix[,r]
            if (use.fast.enrichment.routine == FALSE) {
               GSEA.results <- GSEA.EnrichmentScore(gene.list=gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=correl.matrix[, r])   
            } else {
               GSEA.results <- GSEA.EnrichmentScore2(gene.list=gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=correl.matrix[, r])   
            }
            phi[i, r] <- GSEA.results$ES
        }
        if (fraction < 1.0) { # if resampling then compute ES for all observed rankings
            for (r in 1:nperm) {
                obs.gene.list2 <- obs.order.matrix[,r]
                if (use.fast.enrichment.routine == FALSE) {
                   GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r]) 
                } else {
                   GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])
                }
                obs.phi[i, r] <- GSEA.results$ES
            }
        } else { # if no resampling then compute only one column (and fill the others with the same value)
             obs.gene.list2 <- obs.order.matrix[,1]
            if (use.fast.enrichment.routine == FALSE) {
               GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r]) 
            } else {
               GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])
            }
            obs.phi[i, 1] <- GSEA.results$ES
            for (r in 2:nperm) {
               obs.phi[i, r] <- obs.phi[i, 1]
            }
        }
        gc()
     }

   } else if (reshuffling.type == "gene.labels") { # reshuffling gene labels
      for (i in 1:Ng) {
        gene.set <- gs[[i]]
        #gene.set2 <- vector(length=length(gene.set), mode = "numeric")
        gene.set2<-NULL
		gene.set2 <- which(gene.labels %in% gene.set)
        for (r in 1:nperm) {
            reshuffled.gene.labels <- sample(1:rows)
            if (use.fast.enrichment.routine == FALSE) {
               GSEA.results <- GSEA.EnrichmentScore(gene.list=reshuffled.gene.labels, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.val)   
            } else {
               GSEA.results <- GSEA.EnrichmentScore2(gene.list=reshuffled.gene.labels, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.val)   
            }
            phi[i, r] <- GSEA.results$ES
        }
        if (fraction < 1.0) { # if resampling then compute ES for all observed rankings
           for (r in 1:nperm) {
              obs.gene.list2 <- obs.order.matrix[,r]
              if (use.fast.enrichment.routine == FALSE) {
                 GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
              } else {
                 GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
              }
              obs.phi[i, r] <- GSEA.results$ES
           }
        } else { # if no resampling then compute only one column (and fill the others with the same value)
           obs.gene.list2 <- obs.order.matrix[,1]
           if (use.fast.enrichment.routine == FALSE) {
              GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
           } else {
              GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
           }
           obs.phi[i, 1] <- GSEA.results$ES
           for (r in 2:nperm) {
              obs.phi[i, r] <- obs.phi[i, 1]
           }
        }
        gc()
     }
   }
   
    # Rescaling normalization for each gene set null

   print("Computing rescaling normalization for each gene set null...")
   
   for (i in 1:Ng) {
         pos.phi <- NULL
         neg.phi <- NULL
         for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
               pos.phi <- c(pos.phi, phi[i, j]) 
            } else {
               neg.phi <- c(neg.phi, phi[i, j]) 
            }
         }
         pos.m <- mean(pos.phi)
         neg.m <- mean(abs(neg.phi))

#         if (Obs.ES[i] >= 0) {
#            KS.size[i] <- which.min(abs(KS.mean.table - pos.m))
#         } else {
#            KS.size[i] <- which.min(abs(KS.mean.table - neg.m))
#         }

         pos.phi <- pos.phi/pos.m
         neg.phi <- neg.phi/neg.m
         for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
                phi.norm[i, j] <- phi[i, j]/pos.m
            } else {
                phi.norm[i, j] <- phi[i, j]/neg.m
            }
          }
          for (j in 1:nperm) {
             if (obs.phi[i, j] >= 0) {
                obs.phi.norm[i, j] <- obs.phi[i, j]/pos.m
             } else {
                obs.phi.norm[i, j] <- obs.phi[i, j]/neg.m
             }
          }
          if (Obs.ES[i] >= 0) {
             Obs.ES.norm[i] <- Obs.ES[i]/pos.m
          } else {
             Obs.ES.norm[i] <- Obs.ES[i]/neg.m
          }
   }
   
   
    if(which.pvalue=="NOM p-val") {
      p.vals<-Nominal_pvalue(Ng,phi,nperm,Obs.ES)
	   } else if (which.pvalue=="FWER") {
	    p.vals<-FWER(Ng,nperm,phi.norm,Obs.ES.norm) 
	   } else if (which.pvalue=="FDR q-val") {
	    p.vals<-FDR.mean(Ng,Obs.ES.norm,gs.names,phi.norm,obs.phi.norm,nperm) 
		} else {
		p.vals<-FDR.median(Ng,Obs.ES.norm,gs.names,phi.norm,obs.phi.norm,nperm)
		}
	  return(list(Obs.ES.norm=Obs.ES.norm,p.vals=p.vals))
	}	


#-------------------------------------------------------------------------------------------
# call function Bootstrap_NES
#-------------------------------------------------------------------------------------------

	if(bootstrap==FALSE)
	  {
	     print("Producing result tables and plots...")

         Obs.ES <- signif(Obs.ES, digits=5)
	     report <- data.frame(cbind(gs.names,gs.terms ,size.G,Obs.ES))
         names(report) <- c("GS","TERM", "SIZE", "ES")
         report.index <- order((as.numeric(as.vector(report$ES))), decreasing=TRUE)
         report2<- report[report.index,]
		 report3<-report2[1:topgs,]
	     rownames(report3)<-1:topgs
		 ES.score<-Obs.ES
         time.stamp <- format(Sys.time(), "yy%Ymm%mdd%d hh%Hmm%Mss%S")
		 name<-paste(output.directory,"Summary_",time.stamp,".csv",sep="")
         write.csv(report3, file =name, quote=FALSE, row.names=FALSE)
	  } else {
	    Obs.ES.norm <- vector(length = Ng, mode = "numeric")
        B<-Bootstrap_NES(nperm,Ng,reshuffling.type,gs,gene.labels,gs.names,order.matrix,obs.order.matrix,correl.matrix,obs.correl.matrix,gene.list2,use.fast.enrichment.routine,weighted.score.type,fraction,Obs.ES,which.pvalue)
		Obs.ES.norm<-B$Obs.ES.norm
		p.vals<-B$p.vals
		# Produce results report

       print("Producing result tables and plots...")

       Obs.ES <- signif(Obs.ES, digits=5)
       Obs.ES.norm <- signif(Obs.ES.norm, digits=5)
       p.vals <- signif(p.vals, digits=5)
	   report <- data.frame(cbind(gs.names,gs.terms ,size.G,Obs.ES, Obs.ES.norm, p.vals))
       names(report) <- c("GS","TERM", "SIZE", "ES", "NES", "p-value")
	   if(!is.null(pvalueCutoff))
	   {
         report2 <- report[(as.vector(report$"p-value") < pvalueCutoff),]
         report.index <- order((as.numeric(as.vector(report2$NES))), decreasing=TRUE)
         report2 <- report2[report.index,]
		 reoprt3<-report2[1:topgs,]
	     rownames(report3)<-1:topgs
	   } else {
	     report.index <- order((as.numeric(as.vector(report$NES))), decreasing=TRUE)
         report2 <- report[report.index,]
		 report3<-report2[1:topgs,]
	     rownames(report3)<-1:topgs
		 }
	     ES.score<-Obs.ES.norm
         time.stamp <- format(Sys.time(), "yy%Ymm%mdd%d hh%Hmm%Mss%S")
		 name<-paste(output.directory,"Summary_",time.stamp,".csv",sep="")
         #write.csv(report3, file =name, quote=FALSE, row.names=FALSE)
		 write.table(report3, file = name, sep = ",", row.names = FALSE,
            quote=FALSE,qmethod = "double")
	  
	  }
	 
#---------------------------------------------------------------------------------------------
# GFMaps
# function Genomic.Map
# a) GFMap.Conditioning
# b) Flood_Fill
# c) pvalue_scaling
# d) pvalue_scaling_score
# e) Visualization_Map
#---------------------------------------------------------------------------------------------
# a)	GFMap.Conditioning
#----------------------------------------------------------------------------------------
	
GFMap.Conditioning<-function(gs,ES.score,obs.index,obs.val,A,pval.gene,report3,importance,gs.terms)
 {
 
  #-STEP-1---making a table with following attributes--
  # a) probeset id   b) ID  c) ID_pvalue d)probe_pvalue e)Term
   score<-NULL
   Term_ID<-NULL
   Term<-NULL
   pvalues<-NULL
   for(i in 1:length(obs.index))
     {
         val<-which(obs.index==i)
         score[i]<-obs.val[val]
         Temp<-list(which(lapply(lapply(gs,function(x){x %in% rownames(A)[i]}),sum)==1))
        if(length(Temp[[1]])!=0)
          {
           index<-(Temp[[1]][which(max(ES.score[(Temp[[1]])])==ES.score[Temp[[1]]])[1]])
           Term_ID[i]<-names(index)
           Term[i]<-gs.terms[index]
          } else {
           Term_ID[i]<-NA
           Term[i]<-NA
         }
        pvalues[i]<-pval.gene[val]
     }
     Probeset_id<-rownames(A)
     table<-data.frame(Probeset_id,Term_ID,GeneScore=score,Genepvalue=pvalues,Term)

 
 #-STEP-2-- remove probe sets for which we have no annotation
	GOtable<-table[!(apply(table[2],1,is.na)),]
 

 #--STEP-3--Remove categories with no probes associated with them
	ddm<-NULL
	probes<-NULL
	j<-1
    for(i in 1:dim(report3)[1])
      {
       filter<-function(x) { if(x[2]==report3[i,1]) TRUE else FALSE}
       filter_GO<-apply(GOtable,1,filter)
       if(sum(filter_GO))
         { 
         probes[j]<-sum(filter_GO) 
         ddm<-rbind(ddm,report3[i,])
         j<-j+1
         }
      } 
	ddm<-cbind(ddm,probes)
 #----Rearrangment of the categories depending on importance("ES","Size")
	
  
	if((dim(ddm)[1]) > 11){
	ddf<-ddm[1:11,]
	}else {ddf<-ddm}
	GOtable<-GOtable[(GOtable[,2] %in% ddf[,1]),]
	rownames(ddf)<-ddf[,1]
    if(importance=="size")
	{
     ddf<-ddf[order(ddf$probes,decreasing=TRUE),] 
	} 
  
 #--STEP-4--mapping of probesetids with index of GOterms 
	index<-NULL
	for( i in 1:nrow(GOtable))
		{
			index[i]<-which(GOtable[i,2]==rownames(ddf))
  
		}
	GOtable<-cbind(GOtable,index)
	
	
	return(list(GOtable=GOtable,ddf=ddf,table=table))
	}

#---------------------------------------------------------------------------------------------
# b)	Flood_Fill
#----------------------------------------------------------------------------------------		


Flood_Fill<-function(Sorted_list,relevance)
{
   GOtable<-Sorted_list$GOtable
   ddf<-Sorted_list$ddf
  
 
  
 #--STEP-1--Reorganization of genes based on Significant category & pvalue  
	#------------------------------------------------------------------------------
	 G.Score<-matrix(GOtable$GeneScore,ncol=1)
	 filter<-function(x){ ((x-min(G.Score))/(max(G.Score)-min(G.Score)))}
	 Score.norm<-apply(G.Score,1,filter)
	 GOtable$GeneScore<-Score.norm
	 #-----------------------------------------------------------------------------------
	 scale<-floor(runif(14400,1,(dim(GOtable)[1])))
	 sacle_mat<-GOtable[scale,]
	 GOlist<-vector(mode="list",length=nrow(ddf))
	 if(relevance=="score")
	 {
	 	test1<-as.data.frame(sacle_mat[order(as.numeric(sacle_mat[,6])+as.numeric(sacle_mat[,3]),decreasing=TRUE),])
        for(i in 1:nrow(ddf))
		{
		  val<-(test1[,6]==i)
          sub<-test1[val,c(6,3)]
          GOlist[[i]]<-sub
        }
	 
	 } else {
	 test1<-as.data.frame(sacle_mat[order(as.numeric(sacle_mat[,6])+as.numeric(sacle_mat[,4]),decreasing=FALSE),])
	 for(i in 1:nrow(ddf))
		{
		  val<-(test1[,6]==i)
          sub<-test1[val,c(6,4)]
          GOlist[[i]]<-sub
        }
		
		} 
		 
  #--STEP-2--inital seeds for the Visualization Map
	dimension<-(floor(sqrt(nrow(sacle_mat)))-20)
	mat=matrix(nrow=dimension,ncol=dimension)
    seed1<-c(1,floor(dim(mat)[1]/2),floor(dim(mat)[2]/2))
    seed2<-c(2,floor(dim(mat)[1]/2),floor(dim(mat)[2]*1/6))
    seed3<-c(3,floor(dim(mat)[1]/2),floor(dim(mat)[2]*5/6))
    seed4<-c(4,floor(dim(mat)[1]*1/6),floor(dim(mat)[2]/2))
    seed5<-c(5,floor(dim(mat)[1]*5/6),floor(dim(mat)[2]/2))
    seed6<-c(6,floor(dim(mat)[1]/dim(mat)[1]),floor(dim(mat)[2]/4))
    seed7<-c(7,floor(dim(mat)[1]/dim(mat)[1]),floor(dim(mat)[2]*3/4))
	seed8<-c(8,floor(dim(mat)[1]/4),floor(dim(mat)[2]/1))
	seed9<-c(9,floor(dim(mat)[1]*3/4),floor(dim(mat)[2]/1))
	seed10<-c(10,floor(dim(mat)[1]/1),floor(dim(mat)[2]*3/4))
	seed11<-c(11,floor(dim(mat)[1]/1),floor(dim(mat)[1]*1/4))
	seed12<-c(12,floor(dim(mat)[1]*3/4),floor(dim(mat)[2]/dim(mat)[2]))
	seed13<-c(13,floor(dim(mat)[1]*1/4),floor(dim(mat)[2]/dim(mat)[2]))


	seed<-list(seed1,seed2,seed3,seed4,seed5,seed6,seed7,seed8,seed9,seed10,seed11,seed12,seed13)


	table1<-NULL
	for(i in 1:nrow(ddf))
		{
		table1<-rbind(table1,seed[[i]])
		}
	rownames(table1)<-rownames(ddf)

  #--step-3--initalization of map with significant terms
	Q<-matrix(0,nrow=2)
	sublist<-vector(mode="list",length=nrow(ddf))
	for(i in 1:nrow(ddf))
		{
		x<-as.numeric(table1[i,2])
		y<-as.numeric(table1[i,3]) 
		Q<-cbind(c(x,y))
		sublist[[i]]<-Q
     }
	 
	 
	#--STEP-4---mapping of probesetids->MAP(GOindex)  
		for(i in 1:nrow(sacle_mat))
		{
		index<-as.numeric(sacle_mat[i,6])
		x<-as.numeric(table1[index,2])
		y<-as.numeric(table1[index,3]) 
		Q<-sublist[[index]]
		Q<-matrix(Q,nrow=2)
		if(is.na(mat[x,y]))
			{
			#GOlist[[index]]<-matrix(GOlist[[index]],ncol=2)
			list<-GOlist[[index]][1,]
			mat[x,y]<-list[[1]]+list[[2]]
			GOlist[[index]]<-(GOlist[[index]][-c(1),])
			}
		if(is.na(mat[x-1,y]) && length(mat[x-1,y]==1))
			{
      
			if(!length(which(apply(abs(Q-c(x-1,y)),2,sum) == 0)))
				{ 
					Q<-cbind(Q,c(x-1,y))
				}
			}
		if(x+1 <=dimension && is.na(mat[x+1,y]))
			{
			if (!length(which(apply(abs(Q-c(x+1,y)),2,sum) == 0)))
				{
				Q<-cbind(Q,c(x+1,y))
	             } 
			}
		if(is.na(mat[x,y-1]) && length(mat[x,y-1]==1))
		{  
			if(!length(which(apply(abs(Q-c(x,y-1)),2,sum) == 0)))
				{
				Q<-cbind(Q,c(x,y-1))
				} 
		}
		if(y+1 <=dimension && is.na(mat[x,y+1]))
		{
			if (!length(which(apply(abs(Q-c(x,y+1)),2,sum) == 0)))
				{ 
				Q<-cbind(Q,c(x,y+1))
				}
		}
		Q<-Q[,-c(1)] 
		if(!is.na(Q[1]) && !is.na(Q[2]))
		{
			table1[index,2]<-Q[1]
			table1[index,3]<-Q[2]
       }
		sublist[[index]]<-Q
	}

   
	return(list(mat=mat,ddf=ddf,seed=seed))
	
}
#---------------------------------------------------------------------------------------------
# c)	pvalue_scaling
#----------------------------------------------------------------------------------------	

pvalue_scaling<-function(Template)
{

	mat<-Template$mat
	ddf<-Template$ddf
	seed<-Template$seed
	Sig_matrix<-mat
	
	 for(i in 1:nrow(ddf))
	 {
		x<-which(mat >=i & mat <=i+.01)
		Sig_matrix[x]<-i 
		x<-which(mat >i+.01 & mat <=i+.02)
		Sig_matrix[x]<-i+.2
		x<-which(mat >i+.02 & mat <=i+.03)
		Sig_matrix[x]<-i+.4
		x<-which(mat >i+.03 & mat < i+.05)
		Sig_matrix[x]<-i+.6
		x<-which(mat >=i+.05 & mat < i+1)
		Sig_matrix[x]<-i+.8
	 }
	n<-nrow(ddf)*5
	#Sig_matrix[which(is.na(Sig_matrix))]<-Sig_matrix[which(is.na(Sig_matrix))-3]
	if(max(Sig_matrix)!=nrow(ddf)+.8)
	{
		y<-which(max(Sig_matrix)==Sig_matrix)
		Sig_matrix[y[length(y)]]<-(nrow(ddf)+.8)
	}
  if(min(Sig_matrix)!=1.0)
 {
   Sig_matrix[1]<-1.0
 } 
	return(list(Sig_matrix=Sig_matrix,ns=n,seed=seed,ddf=ddf))
}

#---------------------------------------------------------------------------------------------
#	d) pvalue_scaling_score
#----------------------------------------------------------------------------------------	
pvalue_scaling_score<-function(Template)
{

	mat<-Template$mat
	ddf<-Template$ddf
	seed<-Template$seed
	Sig_matrix<-mat
	
   for(i in 1:nrow(ddf))
     {

        x<-which(mat <=i+1 & mat >=i+.8)
        Sig_matrix[x]<-i
		x<-which(mat <i+.8 & mat >=i+.6)
		Sig_matrix[x]<-i+.2
		x<-which(mat <i+.6 & mat >=i+.4)
		Sig_matrix[x]<-i+.4
		x<-which(mat <i+.4 & mat >= i+.2)
		Sig_matrix[x]<-i+.6
		x<-which(mat <i+.2 & mat >i)
		Sig_matrix[x]<-i+.8
}
n<-nrow(ddf)*5
	#Sig_matrix[which(is.na(Sig_matrix))]<-Sig_matrix[which(is.na(Sig_matrix))-3]
	if(max(Sig_matrix)!=nrow(ddf)+.8)
	{
		y<-which(max(Sig_matrix)==Sig_matrix)
		Sig_matrix[y[length(y)]]<-(nrow(ddf)+.8)
	}
  if(min(Sig_matrix)!=1.0)
 {
   Sig_matrix[1]<-1.0
 } 
	return(list(Sig_matrix=Sig_matrix,n=n,seed=seed,ddf=ddf))
}

#---------------------------------------------------------------------------------------------
# e) Visualization_Map
#----------------------------------------------------------------------------------------	

Visualization_Map<-function(Template_modified)
{
   Sig_matrix<-Template_modified$Sig_matrix
   n<-Template_modified$n
   seed<-Template_modified$seed
   ddf<-Template_modified$ddf  

	colormap<-c("#FFFFFF","#F5F5F5","#DCDCDC","#D3D3D3","#C0C0C0","#FA8D00","#D67800",
	"#B96901","#A05B01","#884D00","#48BDF4","#3DA3D3","#3895C0","#3081A7",
	"#287092","#66C303","#56A503","#4B8E02","#427F02","#346500","#F600FF",
	"#D603DE","#BA03C1","#A501AB","#8F0294","#FF0000","#E20101","#BB0202",
	"#A40202","#900202","#DAED02","#BDCD01","#ABBA04","#8F9B02","#727C01",
	"#77778E","#5F5F72","#434358","#363651","#222245","#0C00F5","#0D02D0",
	"#0902A6","#07017F","#050064","#43E7C2","#38E7A8","#35B497","#31A087",
	"#2A8E77","#7300F5","#6202CE","#5501B3","#4B029C","#3B027B")

	#-----plotting of the image
	testcol<-NULL               
	for(i in 0:(nrow(ddf)-1))
	{
	x1<-colormap[1+(5*i)]
	testcol<-c(testcol,x1)
	}
	col.labels<-rev(rownames(ddf))
	
    return(list(Sig_matrix=Sig_matrix,colormap=colormap,n=n,col.labels=col.labels,testcol=testcol))
	
}
	
#---------------------------------------------------------------------------------------------
#	Genomic.Map
#----------------------------------------------------------------------------------------	
Genomic.Map<-function(gs,ES.score,obs.index,obs.val,A,pval.gene,report3,importance,gs.terms,relevance,output.directory)
{

  Sorted_list<-GFMap.Conditioning(gs,ES.score,obs.index,obs.val,A,pval.gene,report3,importance,gs.terms)
  table<-Sorted_list$table
  Template<-Flood_Fill(Sorted_list,relevance)
  #saveRdata<-paste(output.directory,"Template.RData",sep="")
  #save(Template,file=saveRdata)
  if(relevance=="score"){
  Template_modified<-pvalue_scaling_score(Template)
  } else {
  Template_modified<-pvalue_scaling(Template)
  }
  map<-Visualization_Map(Template_modified)
  Sig_matrix=map$Sig_matrix
  colormap=map$colormap
  n=map$n
  col.labels=map$col.labels
  testcol=map$testcol
  time.stamp <- format(Sys.time(), "yy%Ymm%mdd%d hh%Hmm%Mss%S")
  
  if (.Platform$OS.type == "unix") {
     map.filename <- paste(output.directory,"GFMap_",time.stamp,".pdf", sep="", collapse="")
     pdf(file=map.filename, height = 10.5, width = 10.5)
     } else
	 {
     map.filename <- paste(output.directory,"GFMap_",time.stamp,".pdf", sep="", collapse="")
     pdf(file=map.filename, height = 10, width = 10)
           }
     

	par(mar=c(8,3,4,8))
    image(Sig_matrix,col=colormap[1:n])
    color.legend(1.03,.2,1.06,.9,col.labels,testcol,align="rb",gradient="y")
    par(mar=c(5,4,4,2))
	dev.off()

  
  
  return(list(table=table,Template_modified=Template_modified))
}

    GeneMatrix<-NULL  
	GM<-Genomic.Map(gs,ES.score,obs.index,obs.val,A,pval.gene,report3,importance,gs.terms,relevance,output.directory)	
	GeneMatrix<-GM$table
	
	#---------------------------------------------------------------------
	# Conversion of probe set ID to Entrez Gene ID"
	#----------------------------------------------------------------------
	ent<-paste(affylib,"ENTREZID",sep="")
	entrezIds<-eval(parse(text=paste("entrezIds<-mget(as.character(GeneMatrix[,1]),envir=",ent,")")))
	haveEntrezId <-(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
	GeneMatrix<-GeneMatrix[(GeneMatrix[,1] %in% names(haveEntrezId)),]
	GeneMatrix$Probeset_id<-as.character(haveEntrezId)
	names(GeneMatrix)[1]<-"Enterz.ID"

	Template_modified<-GM$Template_modified
	time.stamp <- format(Sys.time(), "yy%Ymm%mdd%d hh%Hmm%Mss%S")
    filename<-paste(output.directory,"GeneMatrix_",time.stamp,".csv",sep="")
    #write.csv(GeneMatrix, file =filename, quote=FALSE, row.names=FALSE)
	write.table(GeneMatrix, file = filename, sep = ",", row.names = FALSE,
            quote=FALSE,qmethod = "double")
	return(list(GeneMatrix=GeneMatrix,Summary=report3))
	}

