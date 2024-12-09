########Note - The code below closely follows the analysis from https://github.com/micmin3/HPV_OPSCC_Analysis, specifically preproc_celltype.R for preprocessing and cna.R for inferring copy number aberrations
#Load resources
setwd("~/rstuff/resources")
source("functions.R") #Found here: https://github.com/micmin3/HPV_OPSCC_Analysis/blob/main/functions.R
#Colour palettes below - optional
pal1<-readRDS("16pal.RDS")
pal2<-readRDS("35pal.RDS")
pal3<-readRDS("jpal.rds")
#Load data
tpm_mat<-#Load TPM matrix
count_mat<-#Load count matrix
metadata<-#Load metadata table
gtf<-#Load hg38 gtf file
#Add info to metadata file
#Count fraction mitochondrial reads
metadata$mito_frac<-pattern.score(count_mat,"MT-C|MT-N|MT-A|MT-R|MT-T")
#Count fraction noncoding reads
noncoding<-gtf%>%filter(gene_type!="protein_coding")%>%pull(gene_name)%>%unique()
metadata$noncoding_frac<-count_mat[rownames(count_mat)%in%noncoding,]%>%colSums()/colSums(count_mat)
#Count number unique coding genes
metadata$gpc<-apply(count_mat[!rownames(count_mat)%in%noncoding,],2,function(x)sum(x>0))
#Create matrix
m1<-metadata%>%filter(mito_frac<0.3&gpc>=1000)%>%
  rownames()%>%
filter.umi(tpm_mat,.,centre = "none",sparse = T,start="tpm",min_umi = 5000,min_avlog = 4.5)
#Set marker positivity
markers<-c("SOX2","SOX10","TH","PHOX2B")
fx<-mclapply(markers,function(x){
  a<-ifelse(count_mat[x,]>0,"pos","neg")
})
fx2<-do.call(cbind.data.frame,fx)
colnames(fx2)<-markers
metadata<-cbind(metadata,fx2[rownames(metadata),])
#Create UMAP
d1<-sweep(m1,1,rowMeans(m1))%>%
  standard.umap(.,n_neighbors = 50,spread=5, min_dist = 0.0001)
#Cluster UMAP
cl<-cluster.coord(d1,minpts=10,assign_all=T)
metadata[names(cl),"Cluster"]<-cl
#Diagnostic plots
x<-#variable to plot
  ggplot(d1, aes(x=V1, y=V2, colour=metadata[rownames(d1),x])) +
    geom_point(size=2,alpha=0.2) +
    guides(colour = guide_legend(title=x,override.aes = list(size=6,alpha=1)))+
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme_classic(base_size=30) +scale_colour_manual(values=c(pal1,pal2,pal3))
#Cluster DEG
deg<-DE.genes(m1[!rownames(m1)%in%noncoding,],cl[colnames(m1)],min_fc = 1,pairwise = F)
#Cluster enrichment
e<-deg%>%filter(Cluster==x)%>%pull(Gene)%>%enrich(.,rownames(m1[!rownames(m1)%in%noncoding,]))
e$plot
#Assign cluster type after inspecting "deg" and adding a Celltype column based on enrichment analysis and expertise
cn<-metadata%>%filter(!is.na(Cluster))%>%pull(Cluster)%>%unique()%>%sort()
for(i in 1:length(cn)){
  metadata[!is.na(metadata$Cluster)&metadata$Cluster==cn[i],"clustertype"]<-
    deg[deg$Cluster==cn[i],"Celltype"]%>%unique()
}
#Subcluster TME and extract DE genes
m2<-metadata%>%filter(clustertype!="Neuroendocrine"&!is.na(Cluster))%>%
  rownames()%>%
  filter.umi(tpm_mat,.,centre = "none",sparse = T,start="tpm",min_umi = 5000,min_avlog = 4.5)
d11<-sweep(m2,1,rowMeans(m2))%>%
  standard.umap(.,n_neighbors = 50,spread=5, min_dist = 0.0001)
cl2<-cluster.coord(d11,assign_all=T)
metadata[names(cl2),"Cluster_TME"]<-cl2
deg2<-DE.genes(m2[!rownames(m2)%in%noncoding,],metadata[colnames(m2),"Cluster_TME"],min_fc = 1)
cn<-metadata%>%filter(!is.na(Cluster_TME))%>%pull(Cluster_TME)%>%unique()%>%sort()
for(i in 1:length(cn)){
  metadata[!is.na(metadata$Cluster_TME)&metadata$Cluster_TME==cn[i],"clustertype_TME"]<-
    deg2[deg2$Cluster==cn[i],"Celltype"]%>%unique()
}
metadata$clustertype2<-ifelse(is.na(metadata$clustertype_TME),metadata$clustertype,metadata$clustertype_TME)

#Dot plot of canonical markers
markers<-#Canonical markers
m1<-metadata%>%filter(mito_frac<0.3&
                        gpc>=1000&!is.na(clustertype2)&!clustertype2%in%c("Doublet","LQ"))%>%
  rownames()%>%
  filter.umi(tpm_mat,.,centre = "none",sparse = T,start="tpm",min_umi = 5000,min_avlog = 4.5,whitelist = markers$Gene)
l3<-unique(markers$Celltype)
m2<-count_mat[,colnames(m1)]
g1<-intersect(markers$Gene,rownames(m2))

ttr<-metadata[colnames(m1),"Celltype"]
names(ttr)<-colnames(m1)
dot.plot(g1,m1,m2,ttr,l3,minfrac = 0.35,minexp = 1)

##CNA assignment for PPGL and Schwann cells
hg38<-#Load hg38 genome, found in biomart

refs<-c("Fibroblast","Endothelial","Adrenal_Cortex")
metadata%>%filter(clustertype2%in%refs)%>%group_by(sample)%>%count(clustertype2)%>%filter(n>=10)

refcells<-lapply(unique(refs),function(x){
  pl2<-metadata%>%filter(clustertype2==x)%>%count(sample)%>%filter(n>=15)%>%pull(sample)
  unlist(mclapply(pl2,function(y){
    metadata %>% filter(clustertype2==x&sample==y) %>% rownames()%>%sample(.,15)
  }))
})

pl<-metadata%>%filter(clustertype2%in%c("Neuroendocrine","Schwann_Cell"))%>%count(sample)%>%filter(n>=20)%>%pull(sample)

f1<-lapply(pl,function(pat){
  print(paste0("Starting CNA for ",pat))
  
  ep<-metadata %>% filter(sample==pat&clustertype2%in%c("Neuroendocrine","Schwann_Cell")) %>%
    rownames()
 
  rn<-c(ep,unlist(refcells))
  
  merged<-tpm_mat[rownames(tpm_mat)%in%hg38$symbol,rn]
  
  #Create uncentred matrix with only epithelial and stromal cells
  m1<-filter.umi(matrix=merged,cells=rn,centre="none",log=TRUE,min_avlog = 4.5,min_cells = 5000,start="tpm")
  #Set epithelial and reference cells
  query<-ep
  rcx<-refcells
  
  #Create CNA matrix
  mat<-calc.cna(query=query,genome=hg38,matrix=m1,ref=rcx,startlog=TRUE,per_chr = TRUE,
                noise=0.15,mav_n=100,scale=NULL,top_genes=NULL,range=c(-3,3))
  
  
  
  gg<-NULL
  acut=0.15
  dc=0.2
  tr=1/5
  
  if(length(query)<50){a2<-1}
  if(length(query)>=50){
    #Split into subclones
    
    a2<-cluster.clone(mat,query,by_cell = FALSE,knn=TRUE,cluster_method = "louvain",louvain_k =15,
                      cell_matrix = "numeric",merge = TRUE,armcut=10,name=pat,merge_method = "dist",
                      expcut = 0.1,sdcut = 5,top_region = tr,clonesize = 10,adregion = 2/3,
                      adgenes = "top",dimred = "uwot",genome=hg38,
                      adcut=acut,dcut=dc)
    
    #Assign cancer status per subclone
    if(length(a2)>1 & length(query)>=50){
      ad<-a2$ad
      epmat<-a2$epmat
      tcl<-sort(a2$tcl)
      dm<-a2$distmat
      
      a3<-reassign.clones(mat,query,ad,tcl,distmat=dm,em=epmat,tscore=NULL,define_region = tr,
                          adcut=acut,adgenes="top",adregion=2/3,order_clones = TRUE,genome=hg38,
                          name=pat,refilt=FALSE,gt=NULL,goldstandard=gg,clonesize = 10,spec=2)
      #
      ad<-a3$ad
      epmat<-a3$epmat
      tcl<-sort(a3$clones)
      thr<-a3$thresholds
      defmat<-a3$defmat
      
    }
  }
  
  #Assign status if there is only one subclone
  if(length(a2)<2 | length(query)<50){
    r2<-lapply(rcx,length)
    rcx<-rcx[r2>10]
    
    ax<-define.cna(mat,query,unlist(rcx),tumour_score=NULL,top_cells = 1/4,top_region = tr,
                   gs=gg,no_ts=TRUE,name=pat,print_ref = TRUE,extra_score = NULL,spec_cut = 2,genome=hg38)
    ax$plots
    defmat<-ax$def_mat
    defmat<-defmat[,!colnames(defmat) %in% c("knn","prob")]
    thr<-ax$thresholds
    dmx<-defmat[defmat$Origin!="Reference",]
    tcl<-ifelse(dmx$cna3=="Cancer",paste0(pat,"_subclone_A"),
                ifelse(dmx$cna3=="Normal",paste0(pat,"_Normal"),
                       paste0(pat,"_Unresolved")))
    names(tcl)<-rownames(dmx)
    tcl<-sort(tcl)
    
    br1<-br.fun(genes=rownames(mat),separate = "arm",genome = hg38)
    br_vec<-br1$br_vec
    labels<-br1$bdf$labels
    rm<-br1$rm
    
    if(length(rm)!=0){
      mat2<-mat[-rm,]}else{mat2<-mat}
    #
    
    ax2<-admat.fun(mat2,br_vec,tcl[tcl==paste0(pat,"_subclone_A")],
                   labels,cut=acut,genes="top",region=2/3)
    ad<-ax2$ad
    epmat<-ax2$epmat
  }
  
  tc2<-tcl

  tq<-defmat[defmat$Origin=="Reference","Origin"]
  names(tq)<-defmat[defmat$Origin=="Reference",]%>%rownames()
  tq<-sample(tq,50)
  tcl<-c(tq,tcl)
  tcl<-ifelse(tcl=="Reference","Reference",tcl)
  
  df_annot<-metadata[names(tcl),]%>%mutate(Celltype=ifelse(!clustertype2%in%c("Neuroendocrine","Schwann_Cell"),"Stroma",
                clustertype2))%>%select(Celltype,SOX2,SOX10,TH,PHOX2B)

  #Plot CNA matrix
  invisible(plot.cna2(matrix=mat, cells=names(tcl),separate="chr",order=F,
                      name=pat,row_split=tcl,genome = hg38,
                      right_annotation=rowAnnotation(df=df_annot,
                                                     col=list(
                                                       SOX2=c("pos"="black","neg"="white"),
                                                       SOX10=c("pos"="black","neg"="white"),
                                                       TH=c("pos"="black","neg"="white"),
                                                       PHOX2B=c("pos"="black","neg"="white"),
                                                       Celltype=c("PPGL"="red","Schwann_Cell"="blue","Stroma"="green")
                                                     )),
                      column_title_gp = gpar(fontsize = 20),column_title_rot=90,
                      row_title_gp = gpar(fontsize = 20),row_title_rot=0,
                      heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                                grid_height=unit(2,"cm"),
                                                title_gp=gpar(fontsize=20),
                                                labels_gp=gpar(fontsize=20))))
  
  out<-list(clones=tcl,epmat=epmat,defmat=defmat,admat=ad,thresholds=thr,mat=mat)
})
#Add clone info to metadata
setwd#(directory where output from above is saved)
l1<-list.files(pattern = "clones.RDS")
f1<-mclapply(l1,function(x){
  readRDS(x)$defmat
})
dd<-do.call(rbind,f1)
metadata[rownames(dd),"CNA"]<-dd$cna3

#Plot CNA metrics by cell type and marker expression
dr<-dd[unlist(refcells),]
dd<-rbind(dd[dd$Origin!="Reference",],dr)
dd$Celltype<-metadata[rownames(dd),"Celltype"]
dd$SOX10<-metadata[rownames(dd),"SOX10"]
dd$SOX2<-metadata[rownames(dd),"SOX2"]
rm<-dd%>%filter(Celltype=="Neuroendocrine"&cna3!="Cancer")%>%rownames()
dd<-dd[!rownames(dd)%in%rm,]%>%mutate(Subset=ifelse(cna3=="Cancer"&SOX10=="neg"&SOX2=="neg"&Celltype=="Neuroendocrine","Neuroendocrine SOX2-/SOX10-",
                                                    ifelse(cna3=="Cancer"&SOX10=="pos"&Celltype=="Neuroendocrine","Neuroendocrine SOX10+",
                                                           ifelse(cna3=="Cancer"&SOX2=="pos"&Celltype=="Neuroendocrine","Neuroendocrine SOX2+",
                                                                  ifelse(Celltype=="Schwann_Cell","Schwann_Cell","Stroma")))))
dd%>%select(cnascore,cnacor,Subset)%>%rename('CNA signal'=cnascore,'CNA correlation'=cnacor)%>%melt()%>%
  ggplot(.,aes(x=variable,y=value,fill=Subset))+
  theme_classic(base_size = 20)+geom_boxplot()+scale_fill_manual(values=pal1)+xlab("")

#CNA summary figure
nsplit=5
pl<-metadata%>%filter(CNA=="Cancer"&clustertype2=="Neuroendocrine")%>%pull(sample)%>%unique()
f1<-mclapply(pl,function(x){
  #Take CNA matrix of cancer cells only and mean/gene
  r1<-metadata%>%filter(CNA=="Cancer"&sample==x&Celltype=="Neuroendocrine")%>%rownames()
  m1<-readRDS(paste0("/path/to/cna/outputs",x,"_clones.RDS"))$mat[,r1]%>%rowMeans()
  #Split into chromosomes
  br1<-br.fun(genes=names(m1),separate = "chr",genecut = 50,genome=hg38)
  br_vec<-br1$br_vec
  labels<-br1$bdf$labels
  rm<-br1$rm
  
  if(length(rm)!=0){
    m1<-m1[-rm]}else{m1<-m1}
  names(br_vec)<-names(m1)
  bn<-unique(br_vec)
  #Split each chromosome into 5 parts
  bf<-lapply(bn,function(x){
    a<-ntile(br_vec[br_vec==x],nsplit)
    names(a)<-names(br_vec[br_vec==x])
    a})
  names(bf)<-bn
  #Select chromosomes with missing data to set zero
  missing<-setdiff(seq(1,23),bn)
  #Mean CNA value of each chromosome part
  b<-unlist(lapply(seq(1,23),function(y){
    if(!y%in%missing){
      b1<-m1[br_vec==y]
      b2<-bf[[as.character(y)]]
      bq<-unique(b2)
      unlist(lapply(bq,function(z){
        c1<-b1[names(b2[b2==z])]
        mean(c1)
      }))
    }else{c(0,0,0,0,0)}
  }))
  b
})
ef2<-do.call(cbind,f1)
colnames(ef2)<-pl
rownames(ef2)<-seq(1,nrow(ef2))
re<-rownames(ef2)
rs<-character(length=nrow(ef2))
names(rs)<-re
rc<- split(re, cut(seq_along(re), 23, labels = FALSE))
for(i in 1:length(rc)){
  rs[rc[[i]]]<-LETTERS[i]
}

nam <- rep(" ", nrow(ef2))
nam[seq(1,nrow(ef2), nsplit)] <- seq(1,23)
nam[nam=="23"]<-"X"
rownames(ef2) <- nam

df_mut<-#Load mutations data found in Supplementary table S1.
rownames(df_mut)<-df_mut$sample
df_mut$sample<-NULL
head(df_mut)

ef2<-ef2[,rownames(df_mut)]

rs2<-colnames(ef2)
names(rs2)<-colnames(ef2)

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("steelblue", "white", "darkred"))
Heatmap(t(ef2),cluster_columns = F,cluster_rows =F,col=col_fun,show_row_names = F,
        row_split=factor(rs2,levels=unique(rs2)),
        name = "Inferred Copy Number\nLog-ratio",row_title_rot = 0,
        column_title_gp = gpar(fontsize=20),
        right_annotation =  rowAnnotation(df=df_mut[,c("VHL","chr_1p")],
                                          annotation_name_gp=gpar(fontsize=20), simple_anno_size = unit(1, "cm"),
                                          annotation_legend_param=list(grid_width=unit(1, "cm"),
                                                                       grid_height=unit(1,"cm"),
                                                                       title_gp=gpar(fontsize=20),
                                                                       labels_gp=gpar(fontsize=20)),
                                          col=list(VHL=c("Mut"="hotpink","WT"="darkgreen"),
                                                             chr_1p=c("Lost"="navy","Preserved"="green","Unknown"="grey"))),
        heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                  grid_height=unit(1,"cm"),
                                  title_gp=gpar(fontsize=15),
                                  labels_gp=gpar(fontsize=15)),
        show_column_names = F,column_split = rs,column_title = labels,border=T,
        column_title_side = "bottom",column_gap = unit(0,"mm")
)
dev.off()


###Plot CNA correlation and signal for one illustrative sample
clm<-readRDS("/path/to/clones.RDS output from CNA analysis for one sample")
dm<-clm$defmat
dm$Celltype<-metadata[rownames(dm),"Celltype"]
thr<-clm$thresholds
ggplot(dm,aes(x=cnacor,y=cnascore,colour=Celltype))+geom_point(size=2)+theme_classic(base_size = 20)+
  geom_hline(yintercept = thr[2,2],lty=2,lwd=2)+geom_vline(xintercept = thr[2,1],lty=2,lwd=2)+
  xlab("CNA correlation")+ylab("CNA signal")