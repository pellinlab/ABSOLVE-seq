library(data.table)
library("brglm2")
library("ggeffects")
library("ggplot2")
library("pheatmap")
library("parallel")
library("ggeffects")

# inputDataInfo=data.table(Donor=c("0"),Treat=c("0"),cInd=0,cTot=0)
offrAll=NULL
woutDFAll=NULL
# folders=c("raw", "raw_class1","raw_class1_dedudStringent","raw_class1_dedudLoose", "raw_class1_filtered", "raw_class1_deduplicateS", "raw_class1_filtered_deduplicateS")
# folders=c("raw", "raw_class1","raw_class2", "raw_class3", "raw_class4", "raw_class5", "dedudStringent")
# class_i = 'raw_class1'
# folders=c("raw", class_i, paste(class_i,"_dedudStringent",sep=""),paste(class_i,"_dedudLoose",sep=""), paste(class_i,"_filtered",sep=""), paste(class_i,"_deduplicateS",sep=""), 
# class_i = 'dedudStringent_downSampling'
# folders=c("dedudStringentV3", "dedudStringentV3_downSampling10_rep1", "dedudStringentV3_downSampling10_rep2", "dedudStringentV3_downSampling10_rep3", "dedudStringentV3_downSampling10_rep4", "dedudStringentV3_downSampling10_rep5")

# class_i = 'dedudStringentV3_deduplicateS_downSampling001_effects'
# folders=c("dedudStringentV3_deduplicateS", "dedudStringentV3_deduplicateS_downSampling001_rep1", "dedudStringentV3_deduplicateS_downSampling001_rep2", "dedudStringentV3_deduplicateS_downSampling001_rep3", "dedudStringentV3_deduplicateS_downSampling001_rep4", "dedudStringentV3_deduplicateS_downSampling001_rep5")
# class_i = 'dedudStringentV3_downSampling001_effects'
# folders=c("dedudStringentV3", "dedudStringentV3_downSampling001_rep1", "dedudStringentV3_downSampling001_rep2", "dedudStringentV3_downSampling001_rep3", "dedudStringentV3_downSampling001_rep4", "dedudStringentV3_downSampling001_rep5")

class_i = 'fourStrategies_231214'
folders = c("raw", "raw_deduplicateS", "dedudStringentV3" , "dedudStringentV3_deduplicateS")
for(qq in 1:length(folders)){
    print(folders[qq])
    basedir=file.path("/data/pinello/PROJECTS/2022_04_LVOT/trimmed_fastq_vpooled_all/CRISPRessoBatch_alleles_umi",folders[qq])
    basedirRes=file.path("/data/pinello/PROJECTS/2022_04_LVOT/trimmed_fastq_vpooled_all/CRISPRessoBatch_alleles_umi",folders[qq],"results")
    setwd(basedir)
    files=dir(pattern = "*.txt")
    filesSplit=unique(unlist(lapply(strsplit(files,split = "_"),function(x) paste(x[1:4],sep="_",collapse = "_") )))
    woutDF=NULL               
    for(fileS in 1:length(filesSplit)){
        filesLocus=dir(pattern=filesSplit[fileS])
        dataAll=NULL
        offr=read.csv("/data/pinello/PROJECTS/2022_04_LVOT/dataInfomatrix", header = TRUE)
        for(fileL in 1:length(filesLocus)){
          dataRep = fread(filesLocus[fileL], sep="\t",header = T)
          offr[fileL,"Donor"]=unique(unlist(lapply(strsplit(filesLocus[fileL],split = "_"),function(x) x[7] )))
          offr[fileL,"Treat"]=unique(unlist(lapply(strsplit(filesLocus[fileL],split = "_"),function(x) x[8] )))
          offr$cInd[fileL]=sum(dataRep$'#Reads'[dataRep$Unedited=="FALSE"])
          offr$cTot[fileL]=sum(dataRep$'#Reads')
        }
        # offr[,S:=NULL]
        # addind xs, col for graphical output
        offr$xs=as.numeric(factor(offr$Donor))
        offr$xs[offr$Treat=="HiFi"]=offr$xs[offr$Treat=="HiFi"]+0.1
        offr$xs[offr$Treat=="NoEP"]=offr$xs[offr$Treat=="NoEP"]-0.1
        offr$col=ifelse(offr$Treat=="HiFi","cornflowerblue","grey40")

        # make sure that NoEP is the reference level
        offr$Treat=factor(offr$Treat)
        offr$Treat=relevel(offr$Treat,"NoEP")
        # prepare data for logistic regression using success/failures formulation  
        offr$cNet=offr$cTot-offr$cInd

        # run the logistic regression using Biased reduction formula   
        mme=glm(cbind(cInd, cNet) ~ Treat+ Donor, data = offr, family = binomial,method = "brglmFit")
        paramEstTest=as.vector(coef(summary(mme)))
                                                              
        #  glmer(cbind(cInd,cNet) ~ Treat+ (1|Donor), data = offr, family = binomial)
        # given the regression model estimates, calculate the predicted probabilities for Donor and Treat variables and their variances, conf. intervals 
        
        h=ggpredict(mme,c("Donor","Treat"))
        # print(h)
        for(oo in 1:length(h$predicted)){
          offr$predicted[offr$Treat== h$group[oo] & offr$Donor== h$x[oo]]=h$predicted[oo]
          offr$std.error[offr$Treat== h$group[oo] & offr$Donor== h$x[oo]]=h$std.error[oo]
          offr$conf.high[offr$Treat== h$group[oo] & offr$Donor== h$x[oo]]=h$conf.high[oo]
          offr$conf.low[offr$Treat== h$group[oo] & offr$Donor== h$x[oo]]=h$conf.low[oo]
        }
        # calculate marginal proportion (you can imagine this as weighted average across donor)
        he=ggeffect(mme,c("Treat"))       
        # print('marginal')
        # # print(he)
        # for(oo in 1:length(he$predicted)){
        #   print(oo)
        #   print(he$predicted[oo])
        #   print(he$conf.low[oo])
        #   print(he$conf.high[oo])
        # }
        # print('conf.high')
        # print(he$conf.high[1])
        # print(he$conf.high[2])
        # print('conf.low')
        # print(he$conf.low[1])
        # print(he$conf.low[2])                               
        #####################  
        expCoeff=c(h$predicted,h$std.error,h$conf.high,h$conf.low)
        # baselineWT: how many donors have the Treated upper confidence interval execeding the corresponding value in the untreated condition+0.002 (0.2%)
        baselineWT=sum(h$conf.high[seq(2, length(h$predicted),by=2)]< h$predicted[seq(1, length(h$predicted),by=2)]+0.002)

        # baselineFixed: how many donors have the Treated upper confidence interval execeding 0.002 (0.2%)
        baselineFixed=sum(h$conf.high[seq(2, length(h$predicted),by=2)]<0.002)

        # Treat variable coefficient estimate
        pvalueTreatEst=coef(summary(mme))["TreatHiFi","Estimate"]

        # Treat variable p-value  
        # Unilater test b>0
        pvalueTreatPval=(1-pnorm( coef(summary(mme))["TreatHiFi","z value"]))

        # Is Treat variable significant 
        pvalueSigEffectPos=(pvalueTreatEst<0.05)*(pvalueTreatEst>0)

        # Save results in woutDF data.frame
                                
        wout=c(pvalueTreatEst,pvalueTreatPval,paramEstTest,expCoeff,baselineWT,baselineFixed, pvalueSigEffectPos, he$predicted, he$conf.low, he$conf.high)
        if(is.null(woutDF)){
            woutDF=matrix(0,nrow = length(filesSplit), ncol=length(wout))
        }
        # print(filesSplit[fileS])
        # print(wout)
        woutDF[fileS,]=wout
        ndonor=length(unique(offr$Donor))
        offr$ys=offr$cInd/offr$cTot             
        offr$pvalueTreatPval = pvalueTreatPval
        offr$pvalueTreatEst = pvalueTreatEst                     
        offr$pvalueSigEffectPos = pvalueSigEffectPos
                                                   
#         png(filename = file.path(basedirRes,paste(filesSplit[fileS],
#                                                   "IndelsLogisticPen",
#                                                   paste(ifelse((pvalueTreatPval<0.05),"Sig","")),
#                                                   paste(ifelse((pvalueTreatEst>0),"Pos","")),
#                                                   "bWT",baselineWT,
#                                                   "bF",baselineFixed,
#                                                   "SigPos",pvalueSigEffectPos,
#                                                   ".png",sep = "_")),    width =4, height = 4, units = "in", pointsize = 12,res = 300)
        
#         plot(0,0,xlim=c(0.5,ndonor+0.5),ylim=c(0,max(c(0.003,h$conf.high,offr$ys) )),col="white",xlab="Donor",ylab="Estimate proportion",main=paste(filesSplit[fileS],ifelse((pvalueTreatPval<0.05),"*","")),xaxt="n",
#                  sub =paste("Treat. Eff.:", round(summary(mme)$coeff[2,1],4), "Pvalue:",round(summary(mme)$coeff[2,4],4)) )
#         axis(1, at = seq(1, 3, by = 1))
#         abline(v=1:ndonor,lty=2,col="grey80")
#         abline(h=0.001,lty=2,col="red")
#         abline(h=0.002,lty=1,col="red")
#         for(ai in 1:length(offr$conf.low)){
#               segments(offr$xs[ai],
#                        offr$conf.low[ai],
#                        offr$xs[ai],
#                        offr$conf.high[ai],
#                        col=offr$col[ai],
#                        lwd = 5 )
#         }
#         points(offr$xs,offr$predicted,pch="-",col="black",lwd=1)
#         points(offr$xs,offr$ys,pch=20,col="black",cex=1)
#             # Plot the Donor specific value estimated in the untreated + 0.002
#             #points( xs[seq(2, length(h$predicted),by=2)],h$predicted[seq(1, length(h$predicted),by=2)]+0.002,pch=4,col="red",lwd=2)
#         dev.off()
                                                   
        fwrite(offr,file=file.path(basedirRes,filesSplit[fileS]),sep = "\t",col.names = T,row.names = F)
            #####################
            # Generate x values for plotting
        offr$locus=filesSplit[fileS]
        offr$filteringSetting=folders[qq]
        if(is.null(offrAll)){
              offrAll=offr
            }else{
              offrAll=rbind(offrAll,offr)
        }                                 
    }
    woutDF=as.data.frame(woutDF)
    woutDF$offtarget=filesSplit
    woutDF$adjpvalue=p.adjust(woutDF[,2],method = "fdr")
    fwrite(woutDF,file=file.path(basedir,"summaryResultsMargPredict"),sep = "\t",col.names = T,row.names = F)

    sp<-ggplot() +
        geom_point(data=offrAll[offrAll$Treat=="NoEP",], aes(x=factor(Donor), y=ys),position=position_nudge(x = -0.2, y = 0),size = 1,colour="grey40")+
        geom_segment(
          data = offrAll[offrAll$Treat=="NoEP",],
          mapping = aes(x=factor(Donor),
                        y=conf.low,
                        xend=factor(Donor),
                        yend=conf.high
                        #,position=position_dodge(width = 0.50, preserve = "total")
          ),
          position = position_nudge(x = -0.2, y = 0),
          colour="grey40",linewidth = 1
        )+
        geom_point(data=offrAll[offrAll$Treat=="NoEP",], aes(x=factor(Donor), y=predicted),position=position_nudge(x = -0.2, y = 0),size = 1,colour="grey40",shape=3)+
        geom_point(data=offrAll[offrAll$Treat=="HiFi",],
                       aes(x=factor(Donor), y=ys),
                       position=position_nudge(x = 0.2, y = 0),size = 1,colour="cornflowerblue")+
        geom_point(data=offrAll[offrAll$Treat=="HiFi",],
                       aes(x=factor(Donor), y=predicted),
                       position=position_nudge(x = 0.2, y = 0),size = 1,colour="cornflowerblue",shape=3)+

        geom_segment(
              data = offrAll[offrAll$Treat=="HiFi",],
              mapping = aes(x=factor(Donor),
                            y=conf.low,
                            xend=factor(Donor),
                            yend=conf.high
                            #,position=position_dodge(width = 0.50, preserve = "total")
              ),
              position = position_nudge(x = 0.2, y = 0),
              colour="cornflowerblue",linewidth = 1
            )+
            theme_bw(base_size = 8)+
            theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
                  strip.text.x = element_text(angle = 90,vjust = 0.5),
                  strip.background = element_rect(
                    color="white", fill="white")
        )+
        xlab("")+
        ylab("Prop. Indels")+
        #scale_x_discrete("Donor", breaks=c(1:3), labels=c("Donor1","Donor2","Donor3"))+
        geom_vline(xintercept = 1:3, linetype = 3, colour = "grey80")+
        #  ylim(c(0.0,0.1))+
        facet_grid(.~locus)
        ggsave(sp,filename = file.path(basedir,paste( folders[qq],"All.png",sep = "_",collapse = "_")),width = 36,height = 7,units = "in", dpi = 300)

    woutDF$folder=folders[qq]
    woutDFAll=rbind(woutDFAll,woutDF)
}
                                                   
basedir= paste("/data/pinello/PROJECTS/2022_04_LVOT/trimmed_fastq_vpooled_all/CRISPRessoBatch_alleles_umi/globalResults_", class_i, "_compare",sep="")
print(basedir)                     
dir.create(basedir)
fwrite(woutDFAll,file=file.path(basedir,"summaryResultsAllSettings"),sep = "\t",col.names = T,row.names = F)
woutDFAllBack=woutDFAll
woutDFAllBack$marginalPredCrtl=woutDFAll$V46
woutDFAllBack$marginalPredTrear=woutDFAll$V47
woutDFAllBack$confLowCtrl=woutDFAll$V48
woutDFAllBack$confLowTrear=woutDFAll$V49      
woutDFAllBack$confHighCtrl=woutDFAll$V50
woutDFAllBack$confHighTrear=woutDFAll$V51                                  
fwrite(woutDFAllBack,file=file.path(basedir,"summaryResultsAllSettingsPowerMarginal"),sep = "\t",col.names = T,row.names = F)

woutDFAllBack$marginalPredInc=woutDFAllBack$marginalPredTrear-woutDFAllBack$marginalPredCrtl
# woutDFAllBack$
woutDFAll=woutDFAllBack

udonor=unique(offrAll$Donor)
ulocus=unique(offrAll$locus)
ufilteringSetting= unique(woutDFAll$folder)       
                                                   
dd1=matrix(0,nrow = length(ulocus),ncol = length(udonor))
ddA=NULL

for(ufilt in 1:length(ufilteringSetting)){
  dd1=matrix(0,nrow = length(ulocus),ncol = length(udonor))
  for(udon in 1:length(udonor)){
    for(uloc in 1:length(ulocus)){
      #dd1[uloc,udon]=as.numeric(offrAll[Donor==udonor[udon] & Treat=="NoEP" & locus==ulocus[uloc] & filteringSetting==ufilteringSetting[ufilt],"cTot"][1])
      #dd1[uloc,udon]=as.numeric(offrAll[Donor==udonor[udon] & Treat=="NoEP" & locus==ulocus[uloc] & filteringSetting==ufilteringSetting[ufilt],"conf.high"][1])-as.numeric(offrAll[Donor==udonor[udon] & Treat=="NoEP" & locus==ulocus[uloc] & filteringSetting==ufilteringSetting[ufilt],"conf.low"][1])
      #dd1[uloc,udon]=as.numeric(offrAll[Donor==udonor[udon] & Treat=="NoEP" & locus==ulocus[uloc] & filteringSetting==ufilteringSetting[ufilt],"predicted"][1])
      
       dd1[uloc,udon]=as.numeric(offrAll[offrAll$Donor==udonor[udon] & offrAll$Treat=="HiFi" & offrAll$locus==ulocus[uloc] & offrAll$filteringSetting==ufilteringSetting[ufilt],"predicted"][1])-as.numeric(offrAll[offrAll$Donor==udonor[udon] & offrAll$Treat=="NoEP" & offrAll$locus==ulocus[uloc] & offrAll$filteringSetting==ufilteringSetting[ufilt],"predicted"][1])
       dd1[uloc,udon]=ifelse(dd1[uloc,udon]< (0.1/100),0,1)
    } 
  } 
  if(is.null(ddA)){
    ddA=dd1
  }else{
    ddA=cbind(ddA,dd1)
  }
  
}

df=data.frame(Donor=rep(udonor,length(ufilteringSetting)),Filter=as.vector(sapply(ufilteringSetting,function(x) {rep(x,length(udonor))})))
df$Donor_filter=apply(df,1,function(x) paste(x[1],x[2],sep="_"))
rownames(df)=1:nrow(df)

rownames(ddA)=ulocus
colnames(ddA)=df$Donor_filter

# png(filename = file.path(basedir,"proportionHiFiNoEPDiffDiscrete.png"),width =9, height = 9, units = "in", pointsize = 12,res = 600)
# pheatmap(ddA, cluster_rows=F, show_rownames=T,cluster_cols=F,fontsize_row =6,gaps_col=seq(0,ncol(ddA),by=3))
# dev.off()
# which(rownames(ddA) %in% c("1450_OT_0000_REF","1617_OT_0000_REF","1617_OT_0040_ALT"))

png(filename = file.path(basedir,"proportionHiFiNoEPDiffDiscretenoOntarget.png"),  width =9, height = 9, units = "in", pointsize = 12,res = 600)
pheatmap(ddA[-which(rownames(ddA) %in% c("1450_OT_0000_REF","1617_OT_0000_REF","1617_OT_0040_ALT")),], cluster_rows=F, show_rownames=T,cluster_cols=F, fontsize_row =6,gaps_col=seq(0,ncol(ddA),by=3))
# pheatmap(ddA, cluster_rows=F, show_rownames=T,cluster_cols=F, fontsize_row =6,gaps_col=seq(0,ncol(ddA),by=3))
dev.off()

#which(rownames(ddA) %in% c("1450_OT_0000_REF","1617_OT_0000_REF","1617_OT_0040_ALT"))

png(filename = file.path(basedir,"proportionHiFiNoEPDiffDiscrete.png"),    width =9, height = 9, units = "in", pointsize = 12,res = 600)
pheatmap(ddA, cluster_rows=F, show_rownames=T,cluster_cols=F, fontsize_row =6,gaps_col=seq(0,ncol(ddA),by=3))
dev.off()