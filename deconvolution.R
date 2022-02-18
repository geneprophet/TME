##deconvolution
library(readr)
Immune_Landscape <- read_delim("marker/Immune_Landscape.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
length(unique(Immune_Landscape$Symbol))



#############
##SIC 



