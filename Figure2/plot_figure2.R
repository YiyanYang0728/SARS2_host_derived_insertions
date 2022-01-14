library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggtext)

get_boxplot <- function(host, cell_type, str){
        if(host=="human"){
            ref = "Homo_sapiens"
        }else{
            ref = "Chlorocebus_sabaeus"
        }

        infile = paste0(ref, "_chimeric_read_host_genes.tsv")
        df <- read_tsv(infile)
        df <- df %>% filter(Type1=="cdna") %>% separate(Gene, c("Gene", "Version"))

        cell_expr <- read_tsv(paste0(cell_type, "_expr.tsv"))

        if(cell_type=="Vero_cell"){
            data <- df
            cell_expr <- cell_expr %>% mutate(log2expr_1=log2(RPKM_1), log2expr_2=log2(RPKM_2), log2expr_3=log2(RPKM_3)) #%>% mutate(Expr=rowMeans(select(., SARSCOV2_1, SARSCOV2_2, SARSCOV2_3)))
            join_dt <- data %>% left_join(cell_expr, by=c("Gene"="gene_id")) %>% mutate_each(funs(replace(., which(is.na(.)), 0)))
            realExprData1 <- join_dt$log2expr_1
            realExprData2 <- join_dt$log2expr_2
            realExprData3 <- join_dt$log2expr_3
            s1Expr<-cell_expr$log2expr_1
            s2Expr<-cell_expr$log2expr_2
            s3Expr<-cell_expr$log2expr_3
            sel_row = c(1,10,15)
            pal <- c("#DCEEF3", "#005CAB", "#DCEEF3", "#005CAB", "#DCEEF3", "#005CAB")
            plot_dt <- data.frame(Expression_level=c(s1Expr, realExprData1, s2Expr, realExprData2, s3Expr, realExprData3), 
                            Type=c(rep("Background gene\nGSM4916368", length(s1Expr)), 
                                   rep("Gene observed in\nchimeric read\nGSM4916368", length(realExprData1)),
                                   rep("Background gene\nGSM4916369", length(s2Expr)), 
                                   rep("HGene observed in\nchimeric read\nGSM4916369", length(realExprData2)),
                                   rep("Background gene\nGSM4916370", length(s3Expr)), 
                                   rep("Gene observed in\nchimeric read\nGSM4916370", length(realExprData3)))
                            )

        }

        if(cell_type=="Caco_cell"){
            data <- df %>% filter(Cellline=="Caco_cell")
            cell_expr <- cell_expr %>% separate(gene_id, c("gene_id", "version")) %>% mutate(log2expr_1=log2(RPKM_1), log2expr_2=log2(RPKM_2))
            join_dt <- data %>% left_join(cell_expr, by=c("Gene"="gene_id")) %>% mutate_each(funs(replace(., which(is.na(.)), 0)))
            realExprData1 <- join_dt$log2expr_1
            realExprData2 <- join_dt$log2expr_2
            s1Expr<-cell_expr$log2expr_1
            s2Expr<-cell_expr$log2expr_2
            sel_row = c(1,6)
            pal <- c("#DCEEF3", "#005CAB", "#DCEEF3", "#005CAB")
            plot_dt <- data.frame(Expression_level=c(s1Expr, realExprData1, s2Expr, realExprData2), 
                            Type=c(rep("Background gene\nGSM4477888", length(s1Expr)), 
                                   rep("Gene observed in\nchimeric read\nGSM4477888", length(realExprData1)),
                                   rep("Background gene\nGSM4477889", length(s2Expr)), 
                                   rep("Gene observed in\nchimeric read\nGSM4477889", length(realExprData2)))
                            )
        }

        if(cell_type=="Calu_cell"){
            data <- df %>% filter(Cellline=="Calu_cell")
            cell_expr <- cell_expr %>% separate(gene_id, c("gene_id", "version")) %>% mutate(log2expr_1=log2(RPKM_1), log2expr_2=log2(RPKM_2))
            join_dt <- data %>% left_join(cell_expr, by=c("Gene"="gene_id")) %>% mutate_each(funs(replace(., which(is.na(.)), 0)))
            realExprData1 <- join_dt$log2expr_1
            realExprData2 <- join_dt$log2expr_2
            s1Expr<-cell_expr$log2expr_1
            s2Expr<-cell_expr$log2expr_2
            sel_row = c(1,6)
            pal <- c("#DCEEF3", "#005CAB", "#DCEEF3", "#005CAB")
            plot_dt <- data.frame(Expression_level=c(s1Expr, realExprData1, s2Expr, realExprData2), 
                            Type=c(rep("Background gene\nGSM4477962", length(s1Expr)), 
                                   rep("Gene observed in\nchimeric read\nGSM4477962", length(realExprData1)),
                                   rep("Background gene\nGSM4477963", length(s2Expr)), 
                                   rep("Gene observed in\nchimeric read\nGSM4477963", length(realExprData2)))
                            )
        }

        stat.test <- compare_means(
        Expression_level ~ Type, data = plot_dt,
        method = "wilcox.test", alternative = "greater") %>% slice(., sel_row)

        p <- ggboxplot(plot_dt, x="Type", y="Expression_level", fill="Type", palette = pal, legend = "none") +
                    stat_pvalue_manual(stat.test, y.position=17, step.increase=0.06, label="{p.signif}") +
                    labs(title=str, y = "Gene expression level log2(RPKM)", x="") +
                    theme(plot.title = ggtext::element_markdown(), text = element_text(size = 15), axis.text.x=element_text(angle=65, hjust=1))
        return(p)
}

options(repr.plot.width=22, repr.plot.height=12)
p1 <- get_boxplot("human", "Caco_cell", "Homo sapiens Caco-2 cell line")
p2 <- get_boxplot("human", "Calu_cell", "Homo sapiens Calu-3 cell line")
p3 <- get_boxplot("monkey", "Vero_cell", "Chlorocebus sabaeus Vero-6 cell line")

pdf("Figure2.pdf", width=22, height=12)
p1 + p2 + p3 + plot_layout(widths = c(2, 2, 3))
dev.off()
