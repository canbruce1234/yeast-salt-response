# sample usage:
#R --vanilla --args  path/61072_xy.chipchip  path/61072_xy2D.chipchip <chipchip2d3.R

args <- commandArgs();

data1 <- read.table(args[4], header=F, na.strings='\N');
filename=args[4];
cat("filename: ", filename, "\n");
out <- file(args[5],"w");

# 2D local normalization in a square of side 35 (=2*a+1 for a = 17)
a=17

# don't want mitochondrial probes:
idx1=data1[[1]]!='chr17'


chr  <- data1[[1]][idx1];
pos  <- data1[[2]][idx1];
  x  <- data1[[3]][idx1];
  y  <- data1[[4]][idx1];
red  <- data1[[5]][idx1];
green<- data1[[6]][idx1];

 maxy = max(y, na.rm=T)
 y_start= 1:maxy
 y_start=ifelse (y_start<a+1,rep(1,maxy),ifelse (y_start>maxy-a,rep(maxy-2*a,maxy),y_start-a))
 y_end=1:maxy
 y_end=ifelse (y_end<a+1, rep(2*a+1,maxy), ifelse (y_end>maxy-a, rep(maxy,maxy),  y_end+a))
 maxx = max(x, na.rm=T)
 x_start = 1:maxx
 x_start=ifelse (x_start<a+1,rep(1,maxx),ifelse (x_start>maxx-a,rep(maxx-2*a,maxx),x_start-a))
 x_end=1:maxx
 x_end=ifelse (x_end<a+1, rep(2*a+1,maxx), ifelse (x_end>maxx-a, rep(maxx,maxx),  x_end+a))



mred= matrix(rep(NA,maxx*maxy), nrow=maxx, ncol=maxy)
mgreen = mred
mdata=mred

medianred= median(red,na.rm=T)
mediangreen=median(green,na.rm=T)
  for (i in 1:length(chr)) {
    mred[x[i], y[i]] = log2(red[i]/medianred)
    mgreen[x[i], y[i]] = log2(green[i]/mediangreen)
  }

mred[mred > 4.5 | mgreen > 4.5 | mred < -4.5 | mgreen < -4.5] <- NA;
mgreen[mred > 4.5 | mgreen > 4.5 | mred < -4.5 | mgreen < -4.5] <- NA;



# the folowing means, for the case a=10, max(y)=1024
# for y=1-10, starty=1 
# for y=11-1014, i.e., between a+1 and max(y)-a, starty = y-10
# for y=1015-1025, starty=1004, i.e., max(y)-a
 
# The following function "smooth" smooths the matrix of red and green values
# by normalizing them to the local median defined by a square of size
# 2*a+1, centered at a given coordinate point (i.e., the number of 
# data points on either side of the center point is 'a').  The red and green
# values are log normal so they are normalized to a mean log value of 0 with
# a standard deviation (estimated from the median absolute deviation function 
# mad) of 1.


  m2red = mred
  m2green = mgreen
  mdata = mred
  m2data = mdata

  for (i in 1:maxx) {
  if ( i/100 == floor(i/100)  ) {cat("2D smoothing row x=",i,"\n")}
    for (j in 1:maxy) {
        if(is.na(mred[i,j])) {next;}
        medred = median(mred[x_start[i]:x_end[i],y_start[j]:y_end[j]],na.rm=T)
        m2red[i,j] = mred[i,j] - medred;
        medgreen = median(mgreen[x_start[i]:x_end[i],y_start[j]:y_end[j]],na.rm=T)
        m2green[i,j] = mgreen[i,j]- medgreen
    }
  }

mdata=mred-mgreen
m2data=m2red-m2green

  for (i in 1:length(chr)){
     red[i]   =   2^m2red[x[i],y[i]];
     green[i] =   2^m2green[x[i],y[i]];
  }

for (i in 1:length(chr)){
     cat(as.character(chr[i]),"\t", pos[i], "\t", x[i], "\t", y[i], "\t", red[i],"\t", green[i],"\n", sep="", file = out)
}

close(out)

