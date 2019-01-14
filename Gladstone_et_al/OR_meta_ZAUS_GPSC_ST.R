#install.packages("devtools")
#install.packages("metafor")
library("metafor")
library(grid)
library(ggplot2)
library("patchwork")
library(stringr)


#Input T5 from Supplementary T1-T21 in csv format
OR_dat <- read.csv("T5-OR_dataset.csv", header = TRUE, sep =",")

data_subset <- "yeareq_u7_ZAUS_noHIVP"

#Countries that have <5 isolates in the exposed group for GPSC sero analysis are excluded
test_filter <- "c5"
 
country <- unique(OR_dat$Country)
serotypes <- unique(OR_dat$In_Silico_Serotype)

######Serotype and GPSC_serotype##############

#OR per serotype using non-serotype as comparitor (non-exposed)
Sero_res <- matrix(data=NA,nrow=0,ncol=12)
#OR per GPSC_serotype combination using non-GPSC_serotype as comparitor
GPSCsero_res <- matrix(data=NA,nrow=0,ncol=12)
pdf(file="meta_logOR_serotype_GPSC_plots.pdf")
for (sero in  serotypes){
  sero_df <- droplevels(subset(OR_dat, In_Silico_Serotype==sero, select = c("Country","Clinical_Manifest","GPSC_sero")))
  GPSCsero <- unique(sero_df$GPSC_sero)
  sero_tab <- table(subset(OR_dat, In_Silico_Serotype==sero, select = c("Country","Clinical_Manifest")))
  other_sero <- table(subset(OR_dat, In_Silico_Serotype!=sero, select = c("Country","Clinical_Manifest")))
  input <- as.data.frame(cbind(sero_tab,other_sero))
  input <- tibble::rownames_to_column(input,"Country")
  rowcount <- cbind(input,rowSums(input[2:3]))
  input <- subset(rowcount, `rowSums(input[2:3])`>=5)[,1:5]
  colnames(input) <- c("Country","carriage_sero","disease_sero","carriage_nonsero","disease_nonsero")
  if (length(input$Country)>0){
    if (sum(input[,2:3])>0){
      es <- escalc(measure = "OR", ai=disease_sero, bi=carriage_sero, ci=disease_nonsero,di=carriage_nonsero, data = input, append = TRUE, add = 0.5, to = "only0")
      rem <- rma(yi, vi, data= es)
      Sero_res <- rbind(Sero_res,c(sero,"All",sero,signif(predict(rem, transf=exp)$pred, digits = 4),signif(predict(rem, transf=exp)$ci.lb, digits = 4),signif(predict(rem, transf=exp)$ci.ub, digits = 4),colSums(input[2:5]),signif(rem$pval, digits = 2),signif(rem$QEp, digits = 4)))
      forest(rem, slab=es$Country, atransf=exp )
      grid.text(sero, .5, .89, gp=gpar(cex=2))
      grid.text("Log Odds Ratio [95% CI]", .8, .8, gp=gpar(cex=1))
      grid.text(paste("Test of heterogeneity: t2=",signif(rem$tau2, digits = 4), "Q=",signif(rem$QE, digits = 4),"p=",signif(rem$QEp, digits = 4)),  .35, .775, gp=gpar(cex=1))
    }
  }
  for (Gsero in GPSCsero){
    Gsero_tab <- table(subset(OR_dat, GPSC_sero==Gsero, select = c("Country","Clinical_Manifest")))
    other_tab <- table(subset(OR_dat, GPSC_sero!=Gsero, select = c("Country","Clinical_Manifest")))
    input <- as.data.frame(cbind(Gsero_tab,other_tab))
    input <- tibble::rownames_to_column(input,"Country")
    rowcount <- cbind(input,rowSums(input[2:3]))
    input <- subset(rowcount, `rowSums(input[2:3])`>=5)[,1:5]
    colnames(input) <- c("Country","carriage_Gsero","disease_Gsero","carriage_nonGsero","disease_nonGsero")
    if (length(input$Country)>0){
      if (sum(input[,2:3])>0){
        es <- escalc(measure = "OR", ai=disease_Gsero, bi=carriage_Gsero, ci=disease_nonGsero,di=carriage_nonGsero, data = input, append = TRUE, add = 0.5, to = "only0")
        rem <- rma(yi, vi, data= es)
        GPSCsero_res <- rbind(GPSCsero_res,c(sero,Gsero,str_split_fixed(Gsero, "_", 2)[,1],signif(predict(rem, transf=exp)$pred, digits = 4),signif(predict(rem, transf=exp)$ci.lb, digits = 4),signif(predict(rem, transf=exp)$ci.ub, digits = 4),colSums(input[2:5]),signif(rem$pval, digits = 2),signif(rem$QEp, digits = 4)))
        forest(rem, slab=es$Country, atransf=exp )
        grid.text(Gsero, .5, .89, gp=gpar(cex=2))
        grid.text("Log Odds Ratio [95% CI]", .8, .8, gp=gpar(cex=1))
        grid.text(paste("Test of heterogeneity: t2=",signif(rem$tau2, digits = 4), "Q=",signif(rem$QE, digits = 4),"p=",signif(rem$QEp, digits = 4)),  .35, .775, gp=gpar(cex=1))
      }
    }
  }
}
dev.off()

colnames(GPSCsero_res) <- c("Serotype","GPSC_sero", "GPSC","logOR","lower","upper","carriage","disease", "carriage_other", "disease_other","pvalue","heterogeneity")
write.csv(GPSCsero_res, file ="T16_GPSC-serotype_OR.csv", row.names = FALSE)

colnames(Sero_res) <- c("Serotype","GPSC_sero", "GPSC","logOR","lower","upper","carriage","disease", "carriage_other", "disease_other","pvalue","heterogeneity")
write.csv(Sero_res, file ="T15_Serotype_OR.csv", row.names = FALSE)

GPSCsero_res <- rbind(GPSCsero_res,Sero_res)

#Plot OR for each serotype with OR for GPSC_serotype combinations of that serotype
dir.create("per_serotype")
for (sero in serotypes){
  seroplot <- subset(as.data.frame(GPSCsero_res), Serotype==sero, select=c(GPSC,logOR,lower,upper))
  seroplot$logOR <- as.numeric(as.character(seroplot$logOR))
  seroplot$lower <- as.numeric(as.character(seroplot$lower))
  seroplot$upper <- as.numeric(as.character(seroplot$upper))
  if (length(seroplot$logOR) > 2){
    ggplot() +
      ggtitle(sero) +
      geom_pointrange(data=seroplot,mapping=aes(x=GPSC, y=logOR, ymin =lower, ymax =upper),size =0.4) +
      geom_hline(yintercept=1, linetype="dashed", color = "red")  +
      #xlab("GPSC") +
      ylab("OR with 95% CI") +
      #ylab("log OR with 95% CI") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x=element_blank(),
            text = element_text(size=20)) +
      coord_cartesian(ylim=c(0.0, max(seroplot$upper)+5)) +
      #coord_cartesian(ylim=c(0.0,55)) +
      scale_y_continuous(breaks=c(0,1,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,150,200,250,300,350,400,450,500,600,700,800))
    ggsave(paste("per_serotype/",data_subset,"_",sero,"_",test_filter,"_GPSCsero_fittodata.pdf", sep=""), plot= last_plot(), width=5)
    #ggsave(paste("per_serotype/",data_subset,"_",sero,"_",test_filter,"_GPSCsero_standardisedY.pdf", sep=""), plot= last_plot(), width=5)
  }
}

#Fishers comparison of proportion from disease for pairs of GPSCs with n>=5 each of the same serotype
Results_fisher <- matrix(data=NA,nrow=0,ncol=9)
set_sero <- subset(OR_dat, select=c(In_Silico_Serotype,GPSC,Clinical_Manifest,Country))
for (coun in unique(set_sero$Country)){
  coun_tab <- subset(set_sero, Country==coun, select=c(In_Silico_Serotype,Clinical_Manifest,GPSC))
  for (Sero in serotypes){
    sero_tab <- table(subset(coun_tab, In_Silico_Serotype==Sero, select=c(Clinical_Manifest,GPSC)))
    len <- length(colnames(sero_tab))
    if (length(rownames(sero_tab))==1){
      sero_tab <- rbind(sero_tab,0)
      if(!("Disease" %in% row.names(sero_tab))){ 
        row.names(sero_tab) <- c("Carriage","Disease")
      } else {
        row.names(sero_tab) <- c("Disease","Carriage")
      } 
    } else {
    }
    while (length(colnames(sero_tab)) >1) {
      i <- 1
      while (i<length(colnames(sero_tab))){
        i<-i+1
        tab <- rbind(sero_tab[,1],sero_tab[,i])
        if ((sum(tab[1,])>=5) & (sum(tab[2,])>=5)){
          diff <-fisher.test(tab)
          fisher <- diff$p.value
          Results_fisher <- rbind(Results_fisher,c(coun,Sero,colnames(sero_tab)[1],tab[1,"Carriage"],tab[1,"Disease"],colnames(sero_tab)[i],tab[2,"Carriage"],tab[2,"Disease"],fisher))
        }
      }
      sero_tab <- sero_tab[,-1]
    }
  }
}

colnames(Results_fisher) <- c("Country","Sero", "GPSC-A","NP","IPD","GSPC-B","NP","IPD","fisher_p")
write.csv(Results_fisher, file =paste("fishersGPSCpairs",data_subset,"_",test_filter,".csv", sep=""), row.names = FALSE)


#############ST_serotype###################

#OR per ST_serotype combination using non-ST_serotype as comparitor
STsero_res <- matrix(data=NA,nrow=0,ncol=12)
#OR per serotype using non-serotype as comparitor (non-exposed)
pdf(file="meta_logOR_ST_plots.pdf")
for (sero in  serotypes){
  sero_df <- droplevels(subset(OR_dat, In_Silico_Serotype==sero, select = c("Country","Clinical_Manifest","ST_sero")))
  STsero <- unique(sero_df$ST_sero)
    for (Gsero in STsero){
    Gsero_tab <- table(subset(OR_dat, ST_sero==Gsero, select = c("Country","Clinical_Manifest")))
    other_tab <- table(subset(OR_dat, ST_sero!=Gsero, select = c("Country","Clinical_Manifest")))
    input <- as.data.frame(cbind(Gsero_tab,other_tab))
    input <- tibble::rownames_to_column(input,"Country")
    rowcount <- cbind(input,rowSums(input[2:3]))
    input <- subset(rowcount, `rowSums(input[2:3])`>=5)[,1:5]
    colnames(input) <- c("Country","carriage_Gsero","disease_Gsero","carriage_nonGsero","disease_nonGsero")
    if (length(input$Country)>0){
      if (sum(input[,2:3])>0){
        es <- escalc(measure = "OR", ai=disease_Gsero, bi=carriage_Gsero, ci=disease_nonGsero,di=carriage_nonGsero, data = input, append = TRUE, add = 0.5, to = "only0")
        rem <- rma(yi, vi, data= es)
        STsero_res <- rbind(STsero_res,c(sero,Gsero,str_split_fixed(Gsero, "_", 2)[,1],signif(predict(rem, transf=exp)$pred, digits = 4),signif(predict(rem, transf=exp)$ci.lb, digits = 4),signif(predict(rem, transf=exp)$ci.ub, digits = 4),colSums(input[2:5]),signif(rem$pval, digits = 2),signif(rem$QEp, digits = 4)))
        forest(rem, slab=es$Country, atransf=exp )
        grid.text(Gsero, .5, .89, gp=gpar(cex=2))
        grid.text("Log Odds Ratio [95% CI]", .8, .8, gp=gpar(cex=1))
        grid.text(paste("Test of heterogeneity: t2=",signif(rem$tau2, digits = 4), "Q=",signif(rem$QE, digits = 4),"p=",signif(rem$QEp, digits = 4)),  .35, .775, gp=gpar(cex=1))
      }
    }
  }
}
dev.off()

colnames(STsero_res) <- c("Serotype","ST_sero", "ST","logOR","lower","upper","carriage","disease", "carriage_other", "disease_other","pvalue","heterogeneity")
write.csv(STsero_res, file ="T17-ST-Serotype_OR.csv", row.names = FALSE)

STsero_res <- rbind(STsero_res,Sero_res)

#Plot OR for each serotype with OR for ST_serotype combinations of that serotype
for (sero in serotypes){
  seroplot <- subset(as.data.frame(STsero_res), Serotype==sero, select=c(ST,logOR,lower,upper))
  seroplot$logOR <- as.numeric(as.character(seroplot$logOR))
  seroplot$lower <- as.numeric(as.character(seroplot$lower))
  seroplot$upper <- as.numeric(as.character(seroplot$upper))
  if (length(seroplot$logOR) > 2){
    ggplot() +
      ggtitle(sero) +
      geom_pointrange(data=seroplot,mapping=aes(x=ST, y=logOR, ymin =lower, ymax =upper),size =0.4) +
      geom_hline(yintercept=1, linetype="dashed", color = "red")  +
      #xlab("ST") +
      ylab("OR with 95% CI") +
      #ylab("log OR with 95% CI") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x=element_blank(),
            text = element_text(size=20)) +
      coord_cartesian(ylim=c(0.0, max(seroplot$upper)+5)) +
      #coord_cartesian(ylim=c(0.0,55)) +
      scale_y_continuous(breaks=c(0,1,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,150,200,250,300,350,400,450,500,600,700,800))
    ggsave(paste("per_serotype/",data_subset,"_",sero,"_",test_filter,"_STsero_fittodata.pdf", sep=""), plot= last_plot(), width=5)
    #ggsave(paste("per_serotype/",data_subset,"_",sero,"_",test_filter,"_STsero_standardisedY.pdf", sep=""), plot= last_plot(), width=5)
  }
}

Results_fisher <- matrix(data=NA,nrow=0,ncol=9)
set_sero <- subset(OR_dat, select=c(In_Silico_Serotype,ST,Clinical_Manifest,Country))
for (coun in unique(set_sero$Country)){
  coun_tab <- subset(set_sero, Country==coun, select=c(In_Silico_Serotype,Clinical_Manifest,ST))
  for (Sero in serotypes){
    sero_tab <- table(subset(coun_tab, In_Silico_Serotype==Sero, select=c(Clinical_Manifest,ST)))
    len <- length(colnames(sero_tab))
    if (length(rownames(sero_tab))==1){
      sero_tab <- rbind(sero_tab,0)
      if(!("Disease" %in% row.names(sero_tab))){ 
        row.names(sero_tab) <- c("Carriage","Disease")
      } else {
        row.names(sero_tab) <- c("Disease","Carriage")
      } 
    } else {
    }
    while (length(colnames(sero_tab)) >1) {
      i <- 1
      while (i<length(colnames(sero_tab))){
        i<-i+1
        tab <- rbind(sero_tab[,1],sero_tab[,i])
        if ((sum(tab[1,])>=5) & (sum(tab[2,])>=5)){
          diff <-fisher.test(tab)
          fisher <- diff$p.value
          Results_fisher <- rbind(Results_fisher,c(coun,Sero,colnames(sero_tab)[1],tab[1,"Carriage"],tab[1,"Disease"],colnames(sero_tab)[i],tab[2,"Carriage"],tab[2,"Disease"],fisher))
        }
      }
      sero_tab <- sero_tab[,-1]
    }
  }
}

colnames(Results_fisher) <- c("Country","Sero", "ST-A","NP","IPD","GSPC-B","NP","IPD","fisher_p")
write.csv(Results_fisher, file =paste("fishersSTpairs",data_subset,"_",test_filter,".csv", sep=""), row.names = FALSE)
