library(Matrix)
library(ggplot2)
library(ggpubr)
library(mixtools)
library(AdaptGauss)
library(Rtsne)
library(stats)
library(RColorBrewer)
library(dplyr)
library(RSpectra)
library(KRLS)
library(SwarmSVM)

ExtractField<-function(VEC,field,delim){
  return(unlist(strsplit(VEC,split = delim,fixed = T))[field])
}

Read10X<-function (data.dir = NULL, gene.column = 2, unique.features = TRUE) {
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, "barcodes.tsv")
    gene.loc <- file.path(run, "genes.tsv")
    features.loc <- file.path(run, "features.tsv.gz")
    matrix.loc <- file.path(run, "matrix.mtx")
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing")
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing")
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names,
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i],
                                   "_", cell.names)
    }
    feature.names <- read.delim(file = ifelse(test = pre_ver_3,
                                              yes = gene.loc, no = features.loc), header = FALSE,
                                stringsAsFactors = FALSE)
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ",
                    gene.column, " but feature.tsv.gz (or genes.tsv) only has ",
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ",
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[,
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) ==
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls ==
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, ])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data,
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}


Data_Meta<-function(MAT){
  nFeature<-colSums(MAT>0)
  Counts<-colSums(MAT)
  Mt_gene<-grep("^MT-",toupper(rownames(MAT)))
  if(length(Mt_gene)>1){
    Percent.mt<-colSums(MAT[Mt_gene,])/colSums(MAT)
    Meta<-list(nFeature=nFeature,Counts=Counts,Percent.mt=Percent.mt)
  }else if(length(Mt_gene)==1){
    Percent.mt<-MAT[Mt_gene,]/colSums(MAT)
    Meta<-list(nFeature=nFeature,Counts=Counts,Percent.mt=Percent.mt)
  }else{
    Meta<-list(nFeature=nFeature,Counts=Counts)
  }

  return(Meta)
}

Plot_Meta<-function(Meta){
  for (i in 1:length(Meta)) {
    MAT<-data.frame(Value=Meta[[i]],Property=rep(names(Meta)[i],length(Meta[[i]])))
    test<-ggplot(MAT, aes(x=Property, y=Value)) + geom_violin(trim = T,fill='#F8766D')+theme_minimal()+geom_jitter(shape=16, position=position_jitter(0.3))
    assign(paste("p",i,sep = "_"),test)
  }

  if(length(Meta)==2){
    ggarrange(p_1,p_2,nrow = 1)
  }else{
    ggarrange(p_1,p_2,p_3,nrow = 1)
  }
}


Data_subset<-function(MAT,Meta,
                      nFeature.upper=Inf,nFeature.lower=-Inf,
                      Counts.upper=Inf,Counts.lower=-Inf,
                      Percent.mt.upper=Inf,Percent.mt.lower=-Inf){
  if(length(Meta)>2){
    MAT<-MAT[,Meta[[1]]<nFeature.upper&Meta[[1]]>nFeature.lower&Meta[[2]]<Counts.upper&Meta[[2]]>Counts.lower&Meta[[3]]<Percent.mt.upper&Meta[[3]]>Percent.mt.lower]
  }else{
    MAT<-MAT[,Meta[[1]]<nFeature.upper&Meta[[1]]>nFeature.lower&Meta[[2]]<Counts.upper&Meta[[2]]>Counts.lower]
  }

  return(MAT)
}



NormalizeCount<-function(MAT){
  return(sweep(MAT,2,colSums(MAT),FUN = "/")*1e6)
}


MIN_return<-function(x){
  return(min(x[x>0]))
}

Global_Zcut<-function(MAT) {
  VEC<-apply(MAT, 1, MIN_return)
  VEC<-VEC+rnorm(length(VEC),0,0.0001)
  Zcut_univ<-0
  tryCatch({
    MIN_fit = normalmixEM(log(VEC),k = 2)
    INTER<-Intersect2Mixtures(Mean1 = MIN_fit$mu[1],SD1 = MIN_fit$sigma[1],Weight1 = MIN_fit$lambda[1],
                              Mean2 = MIN_fit$mu[2],SD2 = MIN_fit$sigma[2],Weight2 = MIN_fit$lambda[2])
    Zcut_univ<-INTER$CutX
  }, error=function(e){})

  return(exp(Zcut_univ))
}


BIC_LTMG <- function(y, rrr, Zcut) {
  n <- length(y)

  nparams <- nrow(rrr) * 3-1
  w <- rrr[, 1]
  u <- rrr[, 2]
  sig <- rrr[, 3]
  y0 <- y[which(y >= Zcut)]

  cc <- c()
  for (i in 1:nrow(rrr)) {
    c <- dnorm(y0, u[i], sig[i]) * w[i]
    cc <- rbind(cc, c)
  }
  d <- apply(cc, 2, sum)
  e <- sum(log(d))
  f <- nparams * log(n)-e*2
  return (f)
}

BIC_ZIMG <-function(y,rrr,Zcut){
  y<-y[y>Zcut]
  n<-length(y)
  nparams <- nrow(rrr) * 3-1
  w <- rrr[, 1]
  u <- rrr[, 2]
  sig <- rrr[, 3]
  y0 <- y[which(y >= Zcut)]
  cc <- c()
  for (i in 1:nrow(rrr)) {
    c <- dnorm(y0, u[i], sig[i]) * w[i]
    cc <- rbind(cc, c)
  }
  d <- apply(cc, 2, sum)
  e <- sum(log(d))
  f <- nparams * log(n)-e*2
  return (f)
}


Pure_CDF<-function(Vec){
  ### Vec should be sorted ###
  TEMP<-sort(Vec)
  TOTAL<-length(Vec)
  CDF<-rep(0,length = length(TEMP))
  m<-TEMP[1]
  KEEP<-c(1)
  if(length(TEMP)>1){
    for (i in 2:length(TEMP)) {
      if (TEMP[i]==m) {
        KEEP<-c(KEEP,i)
      }else{
        m<-TEMP[i]
        CDF[KEEP]<-(i-1)/TOTAL
        KEEP<-c(i)
      }
    }
  }
  CDF[KEEP]<-1
  return(CDF)
}


KS_LTMG<-function(y,rrr,Zcut){
  y<-sort(y)
  num_c<-nrow(rrr)
  y[which(y<Zcut)]<-Zcut-2
  y0<-y[which(y>=Zcut)]
  p_x<-rep(0,length(y0))

  for(j in 1:num_c){
    p_x<-p_x+pnorm(y0,mean=rrr[j,2],sd=rrr[j,3])*rrr[j,1]
  }

  p_uni_x<-Pure_CDF(y)
  p_uni_x<-p_uni_x[which(y>=Zcut)]
  return(max(abs(p_x-p_uni_x)))
}


KS_ZIMG<-function(y,rrr,Zcut){
  num_c<-nrow(rrr)
  y0<-y[which(y>=Zcut)]
  y0<-sort(y0)
  p_x<-rep(0,length(y0))

  for(j in 1:num_c){
    p_x<-p_x+pnorm(y0,mean=rrr[j,2],sd=rrr[j,3])*rrr[j,1]
  }

  p_uni_x<-Pure_CDF(y0)
  return(max(abs(p_x-p_uni_x)))
}


State_return<-function(x){
  return(order(x,decreasing = T)[1])
}

MINUS<-function(x,y){
  if(x<y){
    return(0)
  }else{
    return(x-y)
  }
}


Fit_LTMG<- function(x, n, q, k, err = 1e-10) {
  q <- max(q, min(x))
  c <- sum(x < q)
  x <- x[which(x >= q)]
  if (length(x) <= k) {
    warning(sprintf("The length of x is %i. Sorry, too little conditions\n", length(x)))
    return(cbind(0, 0, 0))
  }
  mean <- c()
  for (i in 1:k) {
    mean <- c(mean, sort(x)[floor(i * length(x) / (k + 1))])
  }
  mean[1] <- min(x) - 1  # What is those two lines for?
  mean[length(mean)] <- max(x) + 1  # Without them the result of mean[1] is slightly different.
  p <- rep(1 / k, k)
  sd <- rep(sqrt(var(x)), k)
  pdf.x.portion <- matrix(0, length(x), k)

  for (i in 1:n) {
    p0 <- p
    mean0 <- mean
    sd0 <- sd

    pdf.x.all <- t(p0 * vapply(x, function(x) dnorm(x, mean0, sd0), rep(0, k)))
    pdf.x.portion <- pdf.x.all / rowSums(pdf.x.all)
    cdf.q <- pnorm(q, mean0, sd0)
    cdf.q.all <- p0 * cdf.q
    cdf.q.portion <- cdf.q.all / sum(cdf.q.all)
    cdf.q.portion.c <- cdf.q.portion * c
    denom <- colSums(pdf.x.portion) + cdf.q.portion.c
    p <- denom / (nrow(pdf.x.portion) + c)
    im <- dnorm(q, mean0, sd0) / cdf.q * sd0
    im[is.na(im)] <- 0
    mean <- colSums(crossprod(x, pdf.x.portion) + (mean0 - sd0 * im) * cdf.q.portion.c) / denom
    sd <- sqrt((colSums((x - matrix(mean0, ncol = length(mean0), nrow = length(x),
                                    byrow = TRUE)) ^ 2 * pdf.x.portion) + sd0 ^ 2 * (1 - (q - mean0) / sd0 * im) *
                  cdf.q.portion.c) / denom)
    if (!is.na(match(NaN, sd))) {
      break
    }
    if ((mean(abs(p - p0)) <= err) && (mean(abs(mean - mean0)) <= err) &&
        (mean(abs(sd - sd0)) <= err)) {
      break
    }
  }
  return(cbind(p, mean, sd))
}


LTMG<-function(VEC,Zcut_G,k=5){
  y<-log(VEC)
  y<-y+rnorm(length(y),0,0.0001)
  Zcut<-min(log(VEC[VEC>0]))
  if(Zcut<Zcut_G){
    Zcut<-Zcut_G
  }


  if(all(VEC>Zcut_G)){
    rrr<-matrix(c(1,mean(y[y>=Zcut]),sd(y[y>=Zcut])),nrow = 1,ncol = 3)
    MARK<-BIC_ZIMG(y,rrr,Zcut)
    rrr_LTMG<-rrr
    for (K in 2:(k-1)) {
      tryCatch({
        mixmdl<-normalmixEM(y[y>Zcut],K)
        rrr<-cbind(mixmdl$lambda,mixmdl$mu,mixmdl$sigma)
        TEMP<-BIC_ZIMG(y,rrr,Zcut)
        if(TEMP<MARK){
          rrr_LTMG<-rrr
          MARK<-TEMP
        }
      }, error=function(e){})
    }
    rrr_LTMG<-rbind(c(0,-Inf,0.0001),rrr_LTMG)
  }else{
    MARK<-Inf
    rrr_LTMG<-NULL
    for (K in 2:k){
      tryCatch({
        rrr<-Fit_LTMG(y,100,Zcut,K)
        rrr<-matrix(as.numeric(rrr[!is.na(rrr[,2]),]),ncol=3,byrow=F)
        TEMP<-BIC_LTMG(y,rrr,Zcut)
        #print(TEMP)
        if(TEMP<MARK){
          rrr_LTMG<-rrr
          MARK<-TEMP
        }
      }, error=function(e){})
    }
  }

  rrr_LTMG<-rrr_LTMG[order(rrr_LTMG[,2]),]
  rrr_use<-matrix(as.numeric(rrr_LTMG),ncol=3,byrow=F)

  return(rrr_LTMG)
}


find_intersect<-function(aaa){
  rrr<-c()
  for(i in 2:ncol(aaa))
  {
    f <- function(x) dnorm(x, m=aaa[2,i], sd=aaa[3,i]) * aaa[1,i] - dnorm(x, m=aaa[2,1], sd=aaa[3,1]) * aaa[1,1]
    rrr<-c(rrr,uniroot(f, interval=c(aaa[2,i],aaa[2,1]))$root)
  }
  return(rrr)
}



plot_gene<-function(VEC,Data_LTMG,Zcut=-Inf,breaks0=30,Gene=NA){
  if(log(min(VEC[VEC>0]))>Zcut){
    Zcut<-log(min(VEC[VEC>0]))
  }
  
  aaa<-t(Data_LTMG)
  
  
  y<-VEC*0
  y[VEC>Zcut]<-log(VEC[VEC>Zcut])
  y[VEC<=Zcut]<-aaa[2,1]
  
  
  mm<-min(y)
  MM<-max(y)
  diff_c<-MM-mm
  MM<-MM+1
  mm<-mm-1
  x<-seq(mm,MM,by=(MM-mm)/1000)
  h<-hist(y,breaks=breaks0,plot = F)
  nn<-max(h$counts)
  h0<-hist(y[which(y>Zcut)],breaks=breaks0,plot=F)
  
  if(is.na(Gene)){
    plot(h,xlim=c(mm+1,MM-1),ylim=c(0,nn),border="white",col="lightblue")
  }else{
    plot(h,xlim=c(mm+1,MM-1),ylim=c(0,nn),border="white",col="lightblue",main=Gene)
  }
  
  
  n<-length(y)*(h$breaks[2]-h$breaks[1])
  z0<-rep(0,length(x))
  for(i in 1:ncol(aaa))
  {
    z<-dnorm(x,aaa[2,i],aaa[3,i])*aaa[1,i]*n
    abline(v=aaa[2,i],col=c(i+1),lwd=2)
    points(z~x,type="l",col=c(i+1),lwd=5,lty=2)
    z0<-z0+z
    
  }
}
                               
plot_dot<-function(VEC,Data_LTMG,cell_key,Zcut,Gene=NA){
  y<-log(VEC)

  if(min(log(VEC[VEC>0]))>Zcut){
    Zcut<-min(log(VEC[VEC>0]))
  }

  rrr_use<-Data_LTMG
  y[y<=Zcut]<-rrr_use[1,2]

  y_use<-y[y>Zcut]
  y_value<-NULL
  for (k in 1:nrow(rrr_use)) {
    TEMP<-dnorm(y_use,mean = rrr_use[k,2],sd = rrr_use[k,3])*rrr_use[k,1]
    y_value<-rbind(y_value,TEMP)
  }
  y_state<-rep(0,length(y))
  y_state[y>Zcut]<-apply(y_value,2,State_return)-1

  MAT<-data.frame(value=y,state=y_state,Cell=cell_key)

  MAT_use<-MAT %>% group_by(state,Cell) %>% summarise(Expression=mean(value),COUNT=length(Cell))
  MAT_use<-MAT_use %>% group_by(Cell) %>% mutate(Percentage=COUNT/sum(COUNT))
  MAT_use$Percentage<-MAT_use$Percentage*100
  MAT_use$state<-as.factor(MAT_use$state)

  p<-ggplot(MAT_use)+geom_point(aes(x=Cell,y=state,size=Percentage,color=Expression))
  p<-p+scale_color_gradientn(colors = c("#4575B4", "#91BFDB", "#FEE090", "#FC8D59","#D73027"))
  #p<-p+facet_grid(.~Class)
  p<-p+ theme(legend.title = element_text(size = 10),
              legend.text = element_text(size = 12),
              plot.title = element_text(size=16),
              axis.title=element_text(size=16,face="bold"),
              axis.text.y = element_text(size = 10, face = "bold",hjust =0.5),
              axis.text.x = element_text(angle = 30, hjust = 0.5, face = "bold"))
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank())
  p<-p+theme(axis.line.x = element_line(color="black", size = 1),
             axis.line.y = element_line(color="black", size = 1))
  p<-p+scale_x_discrete(name ="Cell Type")+scale_y_discrete(name ="Gene State")

  if(is.na(Gene)){
    print(p)
  }else{
    print(p+ggtitle(Gene))
  }
}




LTMG_MAT<-function(MAT,Zcut_G,Gene_use,k=5){
  
  LTMG_Res<-list(State=matrix(0,nrow = length(Gene_use),ncol = ncol(MAT)),
                 A_MAT=matrix(0,nrow = length(Gene_use),ncol = k),
                 U_MAT=matrix(0,nrow = length(Gene_use),ncol = k),
                 S_MAT=matrix(0.0001,nrow = length(Gene_use),ncol = k))
  
  SEQ<-floor(seq(from = 1,to = length(Gene_use),length.out = 11))
  
  
  for (i in 1:length(Gene_use)) {
    
    if(i %in% SEQ){
      cat(paste0("Progress:",(grep("T",SEQ==i)-1)*10,"%\n" ))
    }  
    
    VEC<-MAT[Gene_use[i],]
    
    y<-log(VEC)
    y<-y+rnorm(length(y),0,0.0001)
    Zcut<-min(log(VEC[VEC>0]))
    if(Zcut<Zcut_G){
      Zcut<-Zcut_G
    }
    
    if(all(VEC>Zcut_G)){
      rrr<-matrix(c(1,mean(y[y>=Zcut]),sd(y[y>=Zcut])),nrow = 1,ncol = 3)
      MARK<-BIC_ZIMG(y,rrr,Zcut)
      rrr_LTMG<-rrr
      for (K in 2:(k-1)) {
        tryCatch({
          mixmdl<-normalmixEM(y[y>Zcut],K)
          rrr<-cbind(mixmdl$lambda,mixmdl$mu,mixmdl$sigma)
          TEMP<-BIC_ZIMG(y,rrr,Zcut)
          if(TEMP<MARK){
            rrr_LTMG<-rrr
            MARK<-TEMP
          }
        }, error=function(e){})
      }
      rrr_LTMG<-rbind(c(0,-Inf,0.0001),rrr_LTMG)
    }else{
      MARK<-Inf
      rrr_LTMG<-NULL
      for (K in 2:k){
        tryCatch({
          rrr<-Fit_LTMG(y,100,Zcut,K)
          rrr<-matrix(as.numeric(rrr[!is.na(rrr[,2]),]),ncol=3,byrow=F)
          TEMP<-BIC_LTMG(y,rrr,Zcut)
          #print(TEMP)
          if(TEMP<MARK){
            rrr_LTMG<-rrr
            MARK<-TEMP
          }
        }, error=function(e){})
      }
    }
    
    if(min(dim(rrr_LTMG))==1){
      y_state<-rep(1,length(y))
    }else{
      rrr_LTMG<-rrr_LTMG[order(rrr_LTMG[,2]),]
      rrr_use<-matrix(as.numeric(rrr_LTMG),ncol=3,byrow=F)
      
      y_use<-y[y>Zcut]
      y_value<-NULL
      for (k in 1:nrow(rrr_use)) {
        TEMP<-dnorm(y_use,mean = rrr_use[k,2],sd = rrr_use[k,3])*rrr_use[k,1]
        y_value<-rbind(y_value,TEMP)
      }
      y_state<-rep(0,length(y))
      y_state[y>Zcut]<-apply(y_value,2,State_return)-1
    }
    
    
    LTMG_Res[[1]][i,]<-y_state
    LTMG_Res[[2]][i,1:nrow(rrr_LTMG)]<-rrr_LTMG[,1]
    LTMG_Res[[3]][i,1:nrow(rrr_LTMG)]<-rrr_LTMG[,2]
    LTMG_Res[[4]][i,1:nrow(rrr_LTMG)]<-rrr_LTMG[,3]
  }
  
  
  rownames(LTMG_Res[[1]])<-Gene_use
  rownames(LTMG_Res[[2]])<-Gene_use
  rownames(LTMG_Res[[3]])<-Gene_use
  rownames(LTMG_Res[[4]])<-Gene_use
  colnames(LTMG_Res[[1]])<-colnames(MAT)
  LTMG_Res[[1]]<-Matrix(LTMG_Res[[1]],sparse = T)
  return(LTMG_Res)
}




LTMG_tsne<-function(File_LTMG,dims=2,perplexity=30,max_iter=5000,partial_pca = F){
  MAT<-t(as.matrix(File_LTMG$State))
  MAT_copy<-NULL
  ## check duplicate
  if(any(duplicated(MAT))){
    CHECK<-grep("T",duplicated(MAT))
    MAT_copy<-MAT
    MAT<-MAT[-CHECK,]
    Position<-rep(0,length(CHECK))

    for (i in 1:length(CHECK)) {
      test<-apply(MAT,1,FUN=all.equal, current=MAT_copy[CHECK[i],])
      Position[i]<-grep("TRUE",test)
    }
  }

  Data_tsne<-Rtsne(MAT,dims=dims,perplexity=perplexity,max_iter=max_iter,partial_pca = partial_pca)

  MAT_tsne<-Data_tsne$Y
  rownames(MAT_tsne)<-rownames(MAT)

  if(is.null(MAT_copy)){
    File_LTMG$tSNE<-MAT_tsne
  }else{
    MAT_temp<-matrix(0,nrow = nrow(MAT_copy),ncol = 2)
    rownames(MAT_temp)<-rownames(MAT_copy)
    MAT_temp[rownames(MAT_tsne),]<-MAT_tsne

    for (i in 1:length(CHECK)) {
      MAT_temp[CHECK[i],]<-MAT_tsne[Position[i],]
    }

    File_LTMG$tSNE<-MAT_temp
  }
  return(File_LTMG)
}



Plot_Cluster<-function(Data_LTMG,Plot_Legend=FALSE,Plot_Label=FALSE){
  if(is.null(Data_LTMG$cluster)){
    cat("Did not find cell cluster information \n")
    cat("Run LTMG_Cluster \n")
    Data_LTMG<-LTMG_Cluster(Data_LTMG)
  }

  Cell_label<-sort(unique(Data_LTMG$cluster))

  cl<-colorRampPalette(brewer.pal(11,"Spectral"))(length(Cell_label))

  plot(x=Data_LTMG$tSNE[names(Data_LTMG$cluster)[Data_LTMG$cluster==Cell_label[1]],1],
       y=Data_LTMG$tSNE[names(Data_LTMG$cluster)[Data_LTMG$cluster==Cell_label[1]],2],
       pch=18,cex=0.7,col=cl[1],
       xlim=c(min(Data_LTMG$tSNE)-2,max(Data_LTMG$tSNE)+2),
       ylim=c(min(Data_LTMG$tSNE)-2,max(Data_LTMG$tSNE)+2),
       frame.plot = F,
       xlab = "t-SNE 1",
       ylab = "t-SNE 2",
       bty='L')
  for (i in 2:length(Cell_label)) {
    points(x=Data_LTMG$tSNE[names(Data_LTMG$cluster)[Data_LTMG$cluster==Cell_label[i]],1],
           y=Data_LTMG$tSNE[names(Data_LTMG$cluster)[Data_LTMG$cluster==Cell_label[i]],2],
           pch=18,cex=0.7,col=cl[i])
  }
  if(Plot_Legend){
    legend("right",legend = Cell_label, fill= cl)
  }
  if(Plot_Label){
    for (i in 1:length(Cell_label)) {
      text(x=mean(Data_LTMG$tSNE[names(Data_LTMG$cluster)[Data_LTMG$cluster==Cell_label[i]],1]),
           y=mean(Data_LTMG$tSNE[names(Data_LTMG$cluster)[Data_LTMG$cluster==Cell_label[i]],2]),
           Cell_label[i],cex=1.5)
    }
  }
}


LTMG_Cluster<-function(File_LTMG,Cut=1e-5,Cluster_size=10){
  if(is.null(File_LTMG$tSNE)){
    cat("Did not project on lower dimension, call tSNE \n")
    File_LTMG<-LTMG_tsne(File_LTMG = File_LTMG)
  }


  if(is.null(File_LTMG$Z)){
    Data_Dis<-as.matrix(gausskernel(File_LTMG$tSNE,sigma = 1))
    diag(Data_Dis)<-1
    S<-Data_Dis
    L<-diag(colSums(S))-S
    ei = eigs_sym(L, 30,which="SM")
    k<-sum(ei$values<Cut)
    Z<-ei$vectors[,(ncol(ei$vectors)-k+1):ncol(ei$vectors)]

    File_LTMG$Z<-Z
  }



  test_clust <- kmeans(File_LTMG$Z, centers=ncol(File_LTMG$Z),iter.max = 500,nstart = 5)
  VEC1<-test_clust$cluster
  names(VEC1)<-colnames(File_LTMG$State)
  VEC1<-VEC1[-which(VEC1 %in% which(table(VEC1)<Cluster_size))]

  test<-clusterSVM(File_LTMG$tSNE[names(VEC1),],VEC1,centers = length(unique(VEC1)),verbose = 0)
  pred=predict(test,File_LTMG$tSNE)
  names(pred$predictions)<-colnames(File_LTMG$tSNE)
  VEC1<-pred$predictions

  NAME<-names(sort(table(VEC1),decreasing = T))
  VEC_temp<-VEC1*0

  for (i in 1:length(NAME)) {
    VEC_temp[VEC1==NAME[i]]<-i
  }
  VEC1<-VEC_temp

  names(VEC1)<-colnames(File_LTMG$State)
  File_LTMG$cluster<-VEC1
  return(File_LTMG)
}





State_count<-function(x,k){
  temp<-table(x)
  VEC<-rep(0,k+1)
  names(VEC)<-0:k
  VEC[names(temp)]<-temp
  return(VEC)
}


State_max<-function(x){
  x<-x[names(x)>0]
  return(c(which.max(x),max(x)))
}

State_test<-function(x){
  x1<-x[1:(length(x)/2)]
  #x1<-x1/sum(x1)
  x2<-x[(length(x)/2+1):length(x)]
  #x2<-x2/sum(x2)
  x<-rbind(x1,x2)
  test<-fisher.test(x)
  return(test$p.value)
}


LTMG_Diff<-function(Data_LTMG,Label,TOP){
  Label_uniq<-sort(unique(Label))
  LTMG_different<-list()
  for (i in 1:length(Label_uniq)) {
    MAT1<-Data_LTMG[[1]][,names(Label)[(Label==Label_uniq[i])]]
    MAT2<-Data_LTMG[[1]][,names(Label)[(Label!=Label_uniq[i])]]
    MAT1<-MAT1[Matrix::rowSums(MAT1)>0,]
    MAT2<-MAT2[Matrix::rowSums(MAT2)>0,]
    ROW<-intersect(rownames(MAT1),rownames(MAT2))
    MAT1<-MAT1[ROW,]
    MAT2<-MAT2[ROW,]
    k<-max(max(MAT1),max(MAT2))
    TEMP1<-t(apply(MAT1, 1, State_count,k=k))
    ROW<-apply(Matrix(TEMP1[,-1]),1,max)
    ROW<-order(ROW,decreasing = T)[1:100]

    TEMP2<-t(apply(MAT2, 1, State_count,k=k))
    TEMP<-cbind(TEMP1[ROW,],TEMP2[ROW,])
    P_value<-apply(TEMP, 1, State_test)

    ORDER<-order(P_value)[1:TOP]

    LTMG_different[[i]]<-cbind(TEMP[ORDER,],P_value[ORDER])
  }
  names(LTMG_different)<-Label_uniq
  File_LTMG$Diff<-LTMG_different
  return(File_LTMG)
}


Plot_State_Heatmap<-function(File_LTMG){

  Cell_ident<-names(File_LTMG$Diff)
  Gene<-NULL
  for (i in 1:length(Cell_ident)) {
    Gene<-c(Gene,rownames(File_LTMG$Diff[[Cell_ident[i]]]))
  }

  Label_uniq<-names(File_LTMG$Diff)
  Label<-factor(File_LTMG$cluster,levels = Label_uniq)


  MAT<-File_LTMG$State[Gene,names(Label)]
  k<-max(MAT)
  test<-apply(MAT,1,function(x) tapply(x,Label,State_count,k=k))

  N<-length(Gene)/length(Label_uniq)

  VEC<-NULL
  for (i in 1:length(Label_uniq)) {
    VEC<-c(VEC,rep(i,N))
  }

  MAT_use<-matrix(0,nrow=length(Gene),ncol = length(Label_uniq))
  ROWNAME<-NULL
  for (i in 1:length(test)) {
    MAT_temp<-matrix(unlist(test[[i]]),ncol = k+1,byrow = T)
    MAT_temp<-sweep(MAT_temp,1,rowSums(MAT_temp),FUN="/")
    MAXSTATE<-which.max(MAT_temp[VEC[i],-1])
    MAT_use[i,]<-MAT_temp[,MAXSTATE+1]
    ROWNAME[i]<-paste(names(test)[i],MAXSTATE,sep = "_")
  }
  rownames(MAT_use)<-ROWNAME
  colnames(MAT_use)<-as.character(Label_uniq)
  MAT_use<-MAT_use[rowSums(MAT_use)>0,]

  Bucket<-melt(MAT_use)
  colnames(Bucket)<-c("Gene_State","Cell_type","Enrichment")
  #Bucket$Cell_type<-factor(as.character(Bucket$Cell_type),levels = Label_uniq)
  Bucket$Gene_State<-factor(Bucket$Gene_State,levels = rev(unique(rownames(MAT_use))))

  p<-ggplot(Bucket,aes(Cell_type,Gene_State))+geom_tile(aes(fill=Enrichment))
  p<-p+scale_fill_gradientn(colors = c("#4575B4", "#91BFDB", "#FEE090", "#FC8D59","#D73027"),limits=c(0,1))
  #p<-p+facet_grid(.~Class)
  p<-p+ theme(legend.title = element_text(size = 10),
              legend.text = element_text(size = 12),
              plot.title = element_text(size=16),
              axis.title=element_text(size=14,face="bold"),
              axis.text.y = element_text(size = 4, face = "bold",hjust =1,),
              axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank())
  p<-p+scale_x_discrete(name ="Cell Type", limits=Label_uniq)+scale_y_discrete(name ="Gene State")
  p
}







