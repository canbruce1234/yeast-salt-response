# Reads three replicates from a Chip-chip experiment performed on a Nimblegen
# yeast genome array. A 7-probe moving window (375 nucleotides on the genome)
# is used to estimate the signal value of the middle probe and the probability 
# that it is significantly different from 0, using the Wilcoxon test.

args <- commandArgs();

data1 <- read.table(args[4], header=TRUE);
data2 <- read.table(args[5], header=TRUE);
data3 <- read.table(args[6], header=TRUE);
out3 <- file(args[7],"w");

o1 <- order(rank(data1[[1]]), data1[[2]]);
o2 <- order(rank(data2[[1]]), data2[[2]]);
o3 <- order(rank(data3[[1]]), data3[[2]]);

chroms <- data1[[1]][o1];
positions <- data1[[2]][o1];

tred1   <- log( data1[[3]][o1]/median(data1[[3]][o1]) ) / log(2);
tgreen1 <- log( data1[[4]][o1]/median(data1[[4]][o1]) ) / log(2);

tred2   <- log( data2[[3]][o2]/median(data2[[3]][o2]) ) / log(2);
tgreen2 <- log( data2[[4]][o2]/median(data2[[4]][o2]) ) / log(2);

tred3   <- log( data3[[3]][o3]/median(data3[[3]][o3]) ) / log(2);
tgreen3 <- log( data3[[4]][o3]/median(data3[[4]][o3]) ) / log(2);

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

         cat(as.character(chr),"\t",pos,"\t",-log(pval)/log 
(10),"\t",signal,"\n", file = out3);
         }
         else {
                 cat(as.character(chr),"\t",pos,"\t",0,"\t",0,"\n",  
file = out3);
         }
}

