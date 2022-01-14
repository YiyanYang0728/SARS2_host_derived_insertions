library(tidyverse)
library(reshape2)
library(patchwork)

##########################################
######### Figure 1 left part #############
##########################################

plot_one_panel <- function(data, geneDfA, kimdata, fill_col, title_str){
    panel = ggplot() + 
    geom_col(data=data, aes(x=Position, y=scaled_count, fill=struct, color=struct)) + 
    geom_rect(data=geneDfA, xmin=0, xmax=29902, ymin=-0.1, ymax=0,
              color = "black",
              size = 0.2) +
    geom_rect(data = geneDfA, aes(xmin = V4, xmax = V5, ymin = -0.1, ymax = 0,
            fill = gene),
            color = "black",
            size = 0.2) +
    geom_segment(data=data, aes(x = 0, y = 0, xend = 29902, yend = 0), arrow = arrow(length = unit(0, "cm"))) +
    scale_fill_manual(values = c(fill_col, "#4373e0", "#979aff", "#eecc00", "#4b4b4b",
                               "#bccee3", "#58bbd5", "#6e25f2", "#ca51f7",
                               "#c72d6e", "#e86501", "#78c200")) + 
    scale_color_manual(values = c(fill_col)) +
    annotate(
        "text", label = "5\'",
        x = -500, y = 0, size = 7
    )+
    annotate(
        "text", label = "3\'",
        x = 30402, y = 0, size = 7
    )+
    scale_y_continuous(breaks=c(0, 0.1974,0.3805,0.5404,0.6747,0.7854,0.8761), limits=c(-1.55, 1.55)) + 
    theme_light(base_size = 18) + 
    # theme(legend.position = "none") +
    ggtitle(title_str)
    return(panel)
}

Structure_Huston<-read.csv("SARS-CoV-2_Full_Length_Secondary_Structure_Map.ct", header=F, sep ="")[-1,]
ifLoop <- function(site){
    flag = "stem"
    if(site==0){return("non-stem")}
    if(site==29903){return("non-stem")}
    cur_pos = Structure_Huston %>% filter(V6==site)
    next_pos = Structure_Huston %>% filter(V6==site+1)
    if((cur_pos$V5==0) || (next_pos$V5==0)){
      flag = "non-stem"
    }
    return(flag)
}

dummy_ct = tibble(Position=1:29903)
dummy_ct = dummy_ct %>% rowwise() %>% mutate(ifloop=ifLoop(Position)) %>% as.data.frame()

# read Human chimeric read jnction sites relative to SARS-CoV-2 reference
dt = read.table("Homo_sapiens_sars2_insert_start.list", col.names=c("Position", "pattern"))
# panel 1
legend = as.character(expression("5'-human-SARS2-3'"))
pattern_ct = dt %>% filter(pattern == "hs") %>% group_by(Position) %>% count() %>% mutate(scaled_count=atan(n/5), chimera_pattern=legend)
ct = left_join(dummy_ct, pattern_ct, by="Position") %>% mutate_if(is.numeric,coalesce,0) %>% 
        mutate(chimera_pattern=legend, struct=ifelse(ifloop=="non-stem","1","0")) %>% arrange(n)
p1 = plot_one_panel(ct, geneDfA, kimdata, c("darkgrey", "#D63729"), paste0("Junction based on ", legend," Chimeric Reads"))

# panel 2
legend = as.character(expression("5'-SARS2-human-3'"))
pattern_ct = dt %>% filter(pattern == "sh") %>% group_by(Position) %>% count() %>% mutate(scaled_count=atan(n/5), chimera_pattern=legend)
ct = left_join(dummy_ct, pattern_ct, by="Position") %>% mutate_if(is.numeric,coalesce,0) %>% 
        mutate(chimera_pattern=legend, struct=ifelse(ifloop=="non-stem","1","0")) %>% arrange(n)
p2 = plot_one_panel(ct, geneDfA, kimdata, c("darkgrey", "#D63729"), paste0("Junction based on ", legend," Chimeric Reads"))
 
# read Monkey chimeric read jnction sites relative to SARS-CoV-2 reference
dt = read.table("Chlorocebus_sabaeus_sars2_insert_start.list", col.names=c("Position", "pattern"))
# panel 3
legend = as.character(expression("5'-monkey-SARS2-3'"))
pattern_ct = dt %>% filter(pattern == "hs") %>% group_by(Position) %>% count() %>% mutate(scaled_count=atan(n/5), chimera_pattern=legend)
ct = left_join(dummy_ct, pattern_ct, by="Position") %>% mutate_if(is.numeric,coalesce,0) %>%
    mutate(chimera_pattern=legend, struct=ifelse(ifloop=="non-stem","1","0")) %>% arrange(n)
p3 = plot_one_panel(ct, geneDfA, kimdata, c("darkgrey", "#D63729"), paste0("Junction based on ", legend," Chimeric Reads"))

# panel 4
legend = as.character(expression("5'-SARS2-monkey-3'"))
pattern_ct = dt %>% filter(pattern == "sh") %>% group_by(Position) %>% count() %>% mutate(scaled_count=atan(n/5), chimera_pattern=legend)
ct = left_join(dummy_ct, pattern_ct, by="Position") %>% mutate_if(is.numeric,coalesce,0) %>% 
    mutate(chimera_pattern=legend, struct=ifelse(ifloop=="non-stem","1","0")) %>% arrange(n)
p4 = plot_one_panel(ct, geneDfA, kimdata, c("darkgrey", "#D63729"), paste0("Junction based on ", legend," Chimeric Reads"))

pdf("Figure_1_left_part.pdf", width=20, height=30)
p1 / p2 / p3 / p4
dev.off()


##########################################
######### Figure 1 right part ############
##########################################

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
    scale_x_continuous(expand = c(1, 1), name ="Number of junction sites in Stems")+
    scale_y_continuous(expand = c(0, 0.1), name = "Percentage (%)")+
    theme_classic()+
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          plot.margin=unit(c(0.3,0.2,0.3,0.2), "cm")) +
    ggtitle(paste(host,"_",mode,"_", pvalueStems, StemsRealData))
    InsertionsInStemsPlot
}

p1 <- plotStemTest("human", "hs")
p2 <- plotStemTest("human", "sh")
p3 <- plotStemTest("monkey", "hs")
p4 <- plotStemTest("monkey", "sh")

pdf("Figure_1_right_part.pdf"), width=15, height=48)
p1 / p2 / p3 / p4
dev.off()
