
## Useful for plotting  via Clarity_Chart
library("Clarity")
## Install using
## remotes::install_github("danjlawson/CLARITY/Clarity")

## Read in the main code
source("mixfunctions.R")
ordermatrix<-function(x,order){
    ## Make a matrix out of x and put the rows and columns in that order
    ## Its really just a coersion to make some functions work properly
    x=x[rownames(x)%in%order,colnames(x)%in%order]
    mat=matrix(NA,nrow=length(order),ncol=length(order))
    rownames(mat)=colnames(mat)=order
    mat[rownames(x),colnames(x)]=x
    mat
}
## Just a heatmap scaling function
myscalefun2=function(x)sign(x)*log(abs(1+x))

################################
################################
################################

## Simulation using two datasets that **share** the mixture edge
set.seed(2)
## Simulate a random graph list with 4 tips, 0 mixtures, and 2 observed datasets
tsim=randomgraphlist(4,0,2,
                     labels=c("A","B","C","O"),outgroup="O")
## Make a similar dataset, but this one containing a mixture edge
tsim1=myregraftmixedgestep(tsim,TRUE)$g
## Extract the parameter vector (NB the parameters are represented on (-inf,inf) and transformed to (0,1) for mixture parameters and (0,inf) for drift parameters.
tsimp1=gparvec(tsim1,TRUE)
np=npars(tsim1)
## Randomize the parameters
tsimp1[1:length(tsimp1)]=rnorm(length(tsimp1))
tsim1[[2]]$mixparmap=1
## Re-parameterise the data
ttruth=parameterise(tsim1,tsimp1,n=2)
## Extract the data matrices
tcov_data=ccov_dag(ttruth)

## Now start from a random graph
tg=randomgraphlist(4,0,2,
                   labels=c("A","B","C","O"),outgroup="O")
## Add a single mixture edge
tgpred=myregraftmixedgestep(tg,TRUE)$g
## Run inference
tg0<- infer_dag(tgpred,tcov_data,maxiter=50,losstol=0.01)
## tg0 is the inferred graph when all parameters are independent

## Duplicate the analysis but fixing the mixture parameters to be the same in each dataset
tg1=tg0$g
## mixparmap tells us which mixture parameters are shared and which are not. The way we've implemented it, we force all graphs to use the parameter of the first graph
tg1[[2]]$mixparmap=1
##Run the inference
tg1<- infer_dag(tg1,tcov_data,maxiter=50,losstol=0.01)
## Extract the predictions
tpred0=ccov_dag(tg0$g)
tpred1=ccov_dag(tg1$g)

## Plotting
options(error = NULL)

par(mfrow=c(3,4))
plot(ttruth,digits=2,main=c("Truth A","Truth B"))
Clarity_Chart(tcov_data[[1]],las=2,text=T,main="Data A")
Clarity_Chart(tcov_data[[2]],las=2,text=T,main="Data B")              
plot(tg0$g,digits=2,main=c("Unpaired A","Unpaired B"))
Clarity_Chart(tpred0[[1]],las=2,text=T,main="Predicted A")
Clarity_Chart(tpred0[[2]],las=2,text=T,main="Predicted B")
plot(tg1$g,digits=2,main=c("Paired A","Paired B"))
Clarity_Chart(tpred1[[1]],las=2,text=T,main="Predicted A")
Clarity_Chart(tpred1[[2]],las=2,text=T,main="Predicted B")

## Numerical comparison of the parameters:
cbind(g1=gparvec(tg1$g),
      g0=gparvec(tg0$g),
      truth=gparvec(ttruth))
## Many problems present!

###############
## Important functions to look at:
## ctree_loss which (in a weird recursive way) computes the sum of the individual graph parameters
## It is just the sum of the squared error between the data and the prediction across all matrices

###############
## Important things to think about:
## What is a good representation of the parameters?
## What is a good representation of the DAG?
## Can we solve for the optima of the parameters more easily? (I believe so; see TreeMix paper)
## Can we rule out entire sub-topologies?
## Can we brute-force examine all graphs?

################################
################################
