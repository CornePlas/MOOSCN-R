## set working directory and define data files
#setwd("/home/corne/Documents/Master-Thesis/R/r/");

## function to check whether required packages are installed
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]); 

## install missing libraries
if (!is.installed('ror')) 
  install.packages("ror", repos="http://cran-mirror.cs.uu.nl/", dependencies=TRUE); 
if (!is.installed('sampling')) 
  install.packages("sampling", repos="http://cran-mirror.cs.uu.nl/", dependencies=TRUE); 
if (!is.installed('chron'))
  install.packages("chron", repos="http://cran-mirror.cs.uu.nl/", dependencies=TRUE); 

## load required libraries
library("ror");
library("sampling");
library("chron");
library("MASS");

file_GA <- "OptionA_nsgaii.txt";

### Read arguments
n <- 5 # if not specified as argument, use 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  n = as.integer(args[1])
}
stopifnot(length(args) < 2)
rm(args)
###

cat(n, file="n.txt")

## read data
data <- as.matrix(read.csv(file_GA, header=FALSE, sep=";"));
data_adj <- 1/data;

## function to find one single solution, it uses random preference statements to filter alternatives
findSol <- function(d, s){
	  pref.s <- c();
	  size <- nrow(d);
	  statements <- 0;
	  iterations <- 0;
	  maxTime <- 0;

	  while( (is.matrix(d)) && (length(d[,1] > 0)) ){ # create loop to end with one variable
	    pref.s <- s * 2; # amount divided by two are the number of statements, thus s
	    selection <- matrix(numeric(1), length(d[,1]), 1);
	    
	    
	    ## determine amount of preference statements
	    x <- 1:length(d[,1]);
	    #x <- x[selection==0];
	    
	    #cat(pref.s, length(d[,1]), pref.s%%2,"\n");
	    while (pref.s > length(d[,1]) || (pref.s == length(d[,1]) && !(pref.s%%2==0))){ # if the remaining alternatives is smaller than the number of alternatives, less statements are needed
	      pref.s <- round((pref.s/2)/2)*2;
	    }
	    cat(format(Sys.time(), "%H:%M;"), " iteration: ",iterations,", size: ",nrow(d),"/",size,", statements:",(pref.s/2)); 
	    cat(format(Sys.time(), "%H:%M;"), " iteration: ",iterations,", size: ",nrow(d),"/",size,", statements:",(pref.s/2), file=paste("output_small-",n,".txt",sep=''), append=TRUE); 
	    statements <- statements + (pref.s / 2); # count no. statements
	    preferences <- matrix(sample(x, pref.s, replace=FALSE), ncol=2, byrow=TRUE);
	    timeIt <- system.time(ror <- utagms(d, preferences, necessary=TRUE, strictVF=TRUE))[1]; # calc necessary relations
	    cat(", time:",timeIt,'\n');
	    cat(", time:",timeIt,'\n',file=paste("output_small-",n,".txt",sep=''),append=TRUE);
	    for (i in 1:length(ror[,1])){ # filter all alternatives that cannot be the preferred one
	      for (j in 1:length(ror[1,])){
		if (i != j && ror[i,j])
		  selection[j] = 1;
	      }
	    }
	    
	    if (timeIt > maxTime){
	      maxTime <- timeIt;
	    }

	    d <- d[selection==0,]; # remove 'dominated' solutions
	    iterations <- iterations + 1;
	  }
	  return(c(d, statements, iterations, maxTime));
}


## use sample
sample_size <- c(10, 20, 30, 40, 50, 60);
no_statements <- c(5,10,20);
results <- matrix(numeric(1), (length(sample_size)*length(no_statements)), 6)
for(j in 1:length(sample_size)){
  for (i in 1:length(no_statements)){
    if (sample_size[j] >= (no_statements[i]*2)){
      sample_data <- srswor(sample_size[j], nrow(data));
      sample_set <- as.matrix(data_adj[sample_data==1,1:3]);
      sys <- system.time(x <- findSol(sample_set, no_statements[i]));
      results[(j-1)*3+i,] <- c(sample_size[j], no_statements[i], sys[1], x[4], x[5], x[6]); 
	# no. alternatives used, time for calc, no. statements, no. iterations, (fixed) statements p/iteration
      cat(sample_size[j],no_statements[i],sys[1],x[4],x[5],x[6],'\n',file=paste("results_small-",n,".txt",sep=''),sep=",",append=TRUE);
    }
  }
}

print(results)
write.matrix(results, file="results_small.csv", sep=";", 1024);

#cat('',file="output_large.txt");
sample_large <- c(80, 100, 120, 140, 160, 180, 200);
statements <- 20;
results <- matrix(numeric(1), 7, 3);
for (i in 1:length(sample_large)){
  sample_data <- srswor(sample_large[i], nrow(data));
  sample_set <- as.matrix(data_adj[sample_data==1,1:3]);
  preferences <- matrix(sample(1:nrow(sample_set), statements, replace=FALSE), ncol=2, byrow=TRUE);
  sys <- system.time(ror <- utagms(sample_set, preferences, necessary=TRUE, strictVF=TRUE))
  results[i,] <- c(sample_large[i], statements, sys[1])
  cat(results[i,],"\n");
  cat("large",results[i,],'\n', file=paste("results_large-",n,".txt",sep=''), sep=",", append=TRUE);
}


print(results)
write.matrix(results, file="results_large.csv", sep=";", 1024);


