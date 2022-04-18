library(tidyverse)
library(reshape2)
library(patchwork)

args = commandArgs(trailingOnly=TRUE)
ref = args[1] # human
pattern = args[2] # hs

bg_struct_annot<-read.table("background.RNA_struct.annot.bed", header=F, sep ="\t")

countStems <- function(Insertions_All){
    count <- 0
    for(i in Insertions_All){
        cur_pos = bg_struct_annot %>% filter(V2==i)
        if((cur_pos$V4 == "S") | (cur_pos$V4 == "S,S")){
          count <- count + 1
        }
    }
    return(count)
}

plotStemTest <- function(host, mode){
    if(host=="human"){
        ref = "Homo_sapiens"
        x_lowlim = 30
        x_uplim = 100
    }else{
        ref = "Chlorocebus_sabaeus"
        x_lowlim = 140
        x_uplim = 400
    }

    infile=paste0(host,".chimeric_junction.RNA_struct.annot.bed")
    dt = read.table(infile, col.names=c("Ref", "Start", "End", "JuncPattern", "Struct"))
    dt <- dt %>% separate(JuncPattern, c("JunctionID", "Pattern"), sep="[.]")

    if(mode=="sh"){
        target_pattern = "sh"
        min_val=1
        max_val=29902-50
    }else{
        target_pattern = "hs"
        min_val=51
        max_val=29902
    }
    dt_1 <- dt %>% filter(Pattern==target_pattern)
    Insertions_All <- dt_1[,2]

    #Real data in Stems
    StemsRealData<-countStems(Insertions_All)
    SimulationStems<-vector()
    for (j in seq(1, 1000)){
    InsPos<-floor(runif(length(Insertions_All), min = min_val, max = max_val))
    SimulationStems<-c(SimulationStems, countStems(InsPos))
    }

    #p-value
    pvalueStems <-length(SimulationStems[SimulationStems < StemsRealData])/1000

    #plot
    InsertionsInStemsPlot<-ggplot()+
    geom_histogram(aes(SimulationStems,
                    after_stat(count*100/sum(count))),
                    fill = "#525252",
                    binwidth = 1, alpha =0.5)+
    geom_vline(xintercept = StemsRealData, color = "#ed1c24", size=3 )+
    scale_x_continuous(expand = c(1, 1), name ="Number of junction sites in stems")+
    scale_y_continuous(expand = c(0, 0.1), name = "Percentage (%)")+
    theme_classic()+
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 30),
          plot.margin=unit(c(0.3,0.2,0.3,0.2), "cm")) +
    ggtitle(paste(host,"_",mode,"_", pvalueStems, StemsRealData))
    InsertionsInStemsPlot
}

p <- plotStemTest(ref, pattern)

pdf(paste0(ref, "_", pattern, ".pdf"), family="ArialMT", width=15, height=12)
p
dev.off()