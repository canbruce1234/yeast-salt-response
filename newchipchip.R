
data1 <- read.table("YAP6/30min/rawdata/47246.clean.new.chipchip", na.strings='\N', header=F);
data2 <- read.table("YAP6/30min/rawdata/57072.clean.new.chipchip", na.strings='\N', header=F);
data3 <- read.table("YAP6/30min/rawdata/47599.clean.new.chipchip", na.strings='\N', header=F);
out3 <- file("YAP6/30min/rawdata/YAP6_30min.clean.exp.out", "w")



o1 <- order(rank(data1[[1]], na.last=NA), data1[[2]]);
o2 <- order(rank(data2[[1]], na.last=NA), data2[[2]]);
o3 <- order(rank(data3[[1]], na.last=NA), data3[[2]]);

chroms <- data1[[1]][o1];
positions <- data1[[2]][o1];

tred1   <- log( data1[[3]][o1]/median(data1[[3]][o1], na.rm=T) ) / log(2);
tgreen1 <- log( data1[[4]][o1]/median(data1[[4]][o1], na.rm=T) ) / log(2);

tred2   <- log( data2[[3]][o2]/median(data2[[3]][o2], na.rm=T) ) / log(2);
tgreen2 <- log( data2[[4]][o2]/median(data2[[4]][o2], na.rm=T) ) / log(2);

tred3   <- log( data3[[3]][o3]/median(data3[[3]][o3], na.rm=T) ) / log(2);
tgreen3 <- log( data3[[4]][o3]/median(data3[[4]][o3], na.rm=T) ) / log(2);


interpol = function(test) {
  if (is.na(test[1])){ test[1] = median(test, na.rm=T)}
  for (i in 2:length(test)) {
    if (is.na(test[i]) & !is.na(test[i+1])) {test[i]=(test[i-1]+test[i+1])/2}
    if (is.na(test[i]) & is.na(test[i+1])) {test[i]=test[i-1]}
  }
test
}


tred1 = interpol(tred1)
tred2 = interpol(tred2)
tred3 = interpol(tred3)
tgreen1 = interpol(tgreen1)
tgreen2 = interpol(tgreen2)
tgreen3 = interpol(tgreen3)

size   <- length(positions);
bandwidth <- 375;
numwidth <- 7;

library(limma);

tcomb <- array(c
(tred1,tgreen1,tred2,tgreen2,tred3,tgreen3),
  dim=c(size,6));
comb <- normalizeQuantiles(tcomb);
red1 <- comb[,1];
green1 <- comb[,2];
red2 <- comb[,3];
green2 <- comb[,4];
red3 <- comb[,5];
green3 <- comb[,6];

for (i in 1:size) {        
         pos <- positions[i];
         chr <- chroms[i];

         if (i >= 2) if (chr != chroms[i-1]) cat (as.character
(chr),"\n");

         for (j in i:max(1,(i-numwidth))) {
                 if (pos - positions[j] < bandwidth & as.character
(chr) == as.character(chroms[j]))
                         left <- j;
         }

         for (j in i:min(size,(i+numwidth))) {
                 if (positions[j] - pos < bandwidth & as.character
(chr) == as.character(chroms[j]))
                         right <- j;
         }

         if (right - left >= 3) {

         red   <- c(red1[left:right], red2[left:right], red3
[left:right]);
         green <- c(green1[left:right], green2[left:right], green3
[left:right]);

         suppressWarnings(ans <- wilcox.test(red, green, paired=T,
alternative = 'greater',conf.int = TRUE));
         pval <- ans$p.value;
         signal <- ans$estimate;

         #cat(pos,"\t",-log(pval)/log(10),"\n", file = out);
         #cat(pos,"\t",signal,"\n", file = out2);
         cat(as.character(chr),"\t",pos,"\t",-log(pval)/log
(10),"\t",signal,"\n", file = out3);
         }
         else {
                 cat(as.character(chr),"\t",pos,"\t",0,"\t",0,"\n",
file = out3);
         }
}

         pos <- positions[i];
         chr <- chroms[i];

         if (i >= 2) if (chr != chroms[i-1]) cat (as.character
(chr),"\n");



############

test=testr1



myfunc = function(test) {
  if (is.na(test[1])){ test[1] = median(test, na.rm=T)}
  for (i in 2:length(test) {
    if (is.na(test[i]) & !is.na(test[i+1])) {test[i]=(test[i-1]+test[i+1])/2}
    if (is.na(test[i]) & is.na(test[i+1])) {test[i]=test[i-1]}
  }
test
}
