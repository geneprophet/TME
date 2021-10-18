#homework2
############
library(readr)
work2 <- read_csv("~/Desktop/助教/work2/work2.csv")
kk = function(x){
  MM <- length(which(x==0))
  MN <- length(which(x==1))
  NN <- length(which(x==2))
  M <- (2*MM + MN)/(2*(MM+MN+NN))
  N <- (2*NN + MN)/(2*(MM+MN+NN))
  EMM <- M*M
  EMN <- 2*M*N
  ENN <- N*N
  p=chisq.test(c(MM, MN, NN), p = c(EMM, EMN, ENN))$p.value
  return(p)
}
work3=work2[,-1]
dt = unlist(apply(work3, 2, kk))
length(which(dt<=0.05))

#homework3
colnames(work2)
df <- work2[,names(which(dt>0.05))]
LD <- cor(df)^2
color <- gray(100:0/100)
library(pheatmap)
pheatmap(LD, cluster_rows = F, cluster_cols = F, col = color, labels_row = "", labels_col = "")


######
df <- read.csv("/Users/kanghongen/Desktop/助教/课件/lesson4/Single.csv")
table(df$rs367896724,df$disease)
chisq.test(table(df$rs367896724,df$disease), correct=F)
X2 = matrix(c(429*2 + 781, 781 + 143*2, 265*2 + 709, 709 + 2 * 177), nrow=2)
X2
chisq.test(X2, correct=F)

df <- read.csv("/Users/kanghongen/Desktop/助教/课件/lesson4/Single.csv")
fit <- glm(disease ~ rs367896724, data=df, family = binomial())
s <- summary(fit)
s$coefficient[2, ]
OR <- exp(s$coefficient[2, 1])
OR

df <- read.csv("/Users/kanghongen/Desktop/助教/课件/lesson4/Multiple.csv")
fit <- glm(disease ~ rs367896724 + rs540431307 + rs555500075, data=df, family = binomial())
s <- summary(fit)
s$coefficient[c(2:4), ]
OR <- exp(s$coefficient[c(2:4), 1])
OR
df <- read.csv("/Users/kanghongen/Desktop/助教/课件/lesson4/Multiple.csv")
fit <- glm(height ~ rs367896724 + rs540431307 + rs555500075, data=df)
s <- summary(fit)
s$coefficient[2, ]

## function for association analysis of continuous traits
getAssoc_c <- function(data, yname, xname, covname = NULL){
  # formula
  if(is.null(covname)){
    fo <- as.formula(paste(yname, "~", xname))
  }else{
    fo <- as.formula(paste(yname, "~", xname, "+", paste(covname, collapse = "+")))
  }
  
  # linear model
  fit <- lm(fo, data = data)
  s <- summary(fit)
  
  # association result
  assocRe <- coef(s)[2, c(1, 2, 4)]
  return(assocRe)
}
## conducting association analysis
# reading data
df <- read.csv("/Users/kanghongen/Desktop/助教/课件/lesson4/Data_for_Association.csv")
snplist <- names(df)[2:101]

# association analysis for continuous traits
yname1 <- "y1"
assocResult1 <- matrix(NA, 100, 4)
assocResult1 <- as.data.frame(assocResult1)
names(assocResult1) <- c("SNP", "Beta", "SE", "P")
Sys.time()
for(i in 1:100){
  xname <- snplist[i]
  assocResult1[i, ] <- c(xname, getAssoc_c(df, yname1, xname))
}
Sys.time()
head(assocResult1)

#########################
dat <-read.table('/Users/kanghongen/Desktop/lesson4/example/extra.ped')
