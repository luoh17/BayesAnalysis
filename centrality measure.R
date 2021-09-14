library("igraph")
library("RColorBrewer")
library("CINNA")
library(StatMatch)
library(igraph)
library(expm)
library(CINNA)
library(markovchain)
library(readxl)
library(quantable)
library(brainGraph)
library(rgl)

randgraph<-function (n,e=n) {
    g2<-random.graph.game(n,e,type="gnm")
    g2<-induced.subgraph(g2,largest_comp(g2))
    g2
  }

largest_comp<-function(graph) {
    cl <- clusters(graph)
    V(graph)[which.max(cl$csize) == cl$membership]
}

getP<-function (g) {
  A<-as.matrix(get.adjacency(g))
  r<-apply(A,1,sum)
  P<-diag(1/r)%*%A
  P
}

pald<-function(X,dist=TRUE,scale=TRUE,show.cont=TRUE,show.plot=TRUE,bet=1,
               gt=NULL,lw=1,tit=FALSE,L=NULL,lab=FALSE,...)
{
  if (dist==FALSE) {
    if (scale==TRUE) { X<-scale(X,scale=TRUE)[,]}
    D<-dist(X) } else D<-X
    D<-as.matrix(D)
    
    getcontmat<-function(D, b=0,h=.5,bet=1,cr=1:dim(D)[1]){
      D<-round(D,15)
      L<-NULL
      n=dim(D)[1]
      A3=matrix(0,n,n)
      for(x in 1:(n-1)){
        for(y in (x+1):n){
          dx=D[x,]; dxt=D[,x]
          dy=D[y,]; dyt=D[,y]
          Uxy=(which((dx<=bet*D[x,y]) | (dy<=bet*D[y,x]))) #the reaching set
          wx<-1*(dx[Uxy]<dy[Uxy])+h*((dx[Uxy]==dy[Uxy]))
          wy<-1*(dy[Uxy]<dx[Uxy])+h*((dx[Uxy]==dy[Uxy]))
          A3[x,Uxy]=A3[x,Uxy]+1/(length(Uxy))*wx
          A3[y,Uxy]=A3[y,Uxy]+1/(length(Uxy))*wy
        }
      }
      diag(A3)<-diag(A3)+b
      rownames(A3)=1:n
      colnames(A3)=1:n
      return(A3/(n-(b==0)))
    }
    
    if (is.null(rownames(D)[1])) {rownames(D)<-1:dim(D)[1]}
    nm<-rownames(D)
    B<-getcontmat(D,h=.5,b=0,bet=bet)
    q<-apply(B,1,sum)
    names(q)<-rownames(D)
    RU<-mean(diag(B))/2
    
    A<-B
    ASym<-pmin(A, t(A))
    
    color<-c( brewer.pal(n = 8, name = "Dark2"), brewer.pal(n=8, name="Set2"))
    color<-rep("grey",16)
    rownames(B)<-rownames(D);colnames(B)<-colnames(D)
    
    diag(ASym)<-1
    g<-graph.adjacency(ASym,weighted=TRUE,mode="undirected")
    g<-simplify(g)
    #V(g)$name<-rownames(D)
    w<-E(g)$weight
    E<-get.edgelist(g)
    E<-E[order(w), ]
    w<-w[order(w)]
    
    
    g<-graph.edgelist(E, directed=FALSE)
    g<-g+setdiff(as.character(1:dim(D)[1]),V(g)$name)
    
    E(g)$weight<-w
    
    Acut<-get.adjacency(g,attr="weight")
    Acut[Acut < RU]<-0
    diag(Acut)<-1
    gcut<-graph.adjacency(Acut,weighted=TRUE,mode="undirected")
    gcut<-simplify(gcut)
    e<-as.numeric(get.edgelist(g)[, 1])
    u<-igraph::clusters(gcut)$membership
    u<-u[order(as.numeric(names(u)))]
    edge_colors<-color[u[e]]
    edge_colors[E(g)$weight<RU]<-"white"
    if(show.cont){
      edge_colors[E(g)$weight<RU]<-"gray"}
    
    
    edge_widths<-E(g)$weight
    edge_widths[edge_widths<RU]<-edge_widths[edge_widths<RU]/20
    
    if (lab) {lab<-rownames(D)[as.numeric(V(g)$name)]} else lab<-""
    
    V(gcut)$clusters<-igraph::clusters(gcut)$membership
    if (is.null(L)) {L<-layout_nicely(g)} else {L<-L[as.numeric(V(g)$name),]}
    if(show.plot){
      plot(g,   ylim=c(-1, 1),xlim=c(-1, 1),...,
           vertex.size=4, vertex.label.cex=1,
           vertex.color=color[igraph::clusters(gcut)$membership],
           #           vertex.label.color=color[igraph::clusters(gcut)$membership],   
           vertex.label.color="black",
           vertex.label=lab,
           vertex.label.dist = 1, edge.width=lw*100*(edge_widths),
           edge.color=edge_colors, asp=0,layout=L,
           main="PaLD")
      if ((!dist)&(tit))
      {title(paste(abbreviate(colnames(X)),collapse=","))}
    }
    if(!is.null(gt[1])){
      plot(g,   ylim=c(-1, 1),xlim=c(-1, 1),...,
           vertex.size=4, vertex.label.cex=1.2,
           #           vertex.color=color[igraph::clusters(gcut)$membership],
           #color[igraph::clusters(gcut)$membership],
           vertex.label.color=gt[as.numeric(V(g)$name)],
           vertex.color=gt[as.numeric(V(g)$name)],
           #vertex.label.color="black",
           vertex.label=lab[as.numeric(V(g)$name)],
           vertex.label.dist = 1, edge.width=lw*100*(edge_widths),
           edge.color=edge_colors, asp=0,layout=L)
      if ((!dist)&(tit))
      {title(paste(abbreviate(colnames(X)),collapse=","))}
      V(g)$gt<-gt[as.numeric(V(g)$name)]
    }
    V(g)$name<-rownames(D)[as.numeric(V(g)$name)]
    
    
    rownames(ASym)<-rownames(B);colnames(ASym)<-colnames(B)
    
    cl<-igraph::clusters(gcut)$membership
    names(cl)<-rownames(D)[as.numeric(names(cl))]
    cl<-cl[sapply(rownames(D),function(c) which(names(cl)==c))]
    
    V(gcut)$name<-rownames(D)[as.numeric(V(gcut)$name)]
    
    list(C=B,Cmin=ASym,g=g,g2=gcut,bound=RU,clusters=(cl),
         isolated=setdiff(rownames(D),V(g)$name),layout=L,depths=(q))

    
    }

### cent_compute ###

exc<-proper_centralities(g1)[c(6,31,33,38,43)]
cent_compute<-function(g1,prcent=setdiff(proper_centralities(g1),exc)){
  # CINNA centrality measures
  print(g1$name)
  M2<-calculate_centralities(g1, include = prcent)
  inc2<-which(as.vector(unlist(lapply(M2,length)))==vcount(g1))
  M<-as.data.frame(M2[inc2])
  colnames(M)<-names(M2[inc2])
  M<-data.frame(M)
  
  P<-getP(g1)
  MC <- new("markovchain",transitionMatrix = P,name = "MC")
  D <- meanFirstPassageTime(MC)
  
  Pald_short<-pald(shortest.paths(g1),dist=TRUE,show.plot=FALSE)$depths # Pald - shortest path
  Pald_hitTo<-pald(D,dist=TRUE,show.plot=FALSE)$depths # Pald - hitting time to
  Pald_hitFrom<-pald(t(D),dist=TRUE,show.plot=FALSE)$depths # Pald - hitting time from
  
  M$Pald_short = Pald_short
  M$Pald_hitTo = Pald_hitTo
  M$Pald_hitFrom = Pald_hitFrom
  
  # rownames(M) <- paste(c(rownames(M),g1$name),collapse=".")
  M <- as.matrix(M)
  M
}


# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}



# ----- END plot function ----- #

####### ------ paldsmooth ------ #######
paldsmooth<-function(M,V,name=rownames(M),diagz=FALSE,R=FALSE,global=TRUE,local=FALSE,original=FALSE){
  M<-cbind(M)
  V<-cbind(V)
  Mt<-M
  D<-gower.dist(Mt)
  L<-pald(D,dist=TRUE,show.plot=FALSE)
  if (diagz) {diag(L$C)<-0;diag(L$Cmin)<-0}
  p<-apply(Mt,1,paste,collapse=",")
  nm<-as.numeric(V(L$g2)$name)
  
  if (R) dev.new()
  
  #  plot(L$g2,vertex.label.color="blue",vertex.label.cex=.8,
  #       layout=L$layout,vertex.size=1,vertex.label=p[nm])
  
  R2<-NULL
  hc<-colorRampPalette(brewer.pal(8, "PiYG"))(40)
  V<-cbind(V)
  for (i in 1:dim(V)[2])
  {
    vc<-V[,i]
    varname<-colnames(V)[i]
    col<-hc[round((vc-min(vc))/(max(vc)-min(vc))*35)+1]
    Cweightg<-diag(1/apply(L$C,1,sum))%*%L$C
    
    vcg<-as.vector(Cweightg%*%vc)  # PaLD smoothing
    colg<-hc[round((vcg-min(vcg))/(max(vcg)-min(vcg))*35)+1]
    Clocal<-L$C
    Clocal[Clocal<L$bound]<-0.001
    Cweightl<-diag(1/apply(Clocal,1,sum))%*%Clocal
    vcl<-as.vector(Cweightl%*%vc)  # PaLD smoothing
    coll<-hc[round((vcl-min(vcl))/(max(vcl)-min(vcl))*35)+1]
    diff<-sum( (vc-mean(vc))^2)
    diffl<-sum( (vc-vcl)^2 )
    diffg<-sum( (vc-vcg)^2 )
    R2<-rbind(R2,c(diff,diffl,diffg,
                   (diff-diffl) /diff,
                   (diff-diffg) /diff))
    titg<-c(round((diff-diffg) /diff,4))
    titl<-c(round((diff-diffl) /diff,4))
    
    if (global)
    {if (R) dev.new()
      plot(L$g2,vertex.label.color="black",vertex.label.cex=0.1,
           layout=L$layout,vertex.size=7,vertex.color=colg[nm],vertex.label="")
      title(c("global",varname,titg))}
    
    
    if (local)
    {if (R) dev.new()
      plot(L$g2,vertex.label.color="black",vertex.label.cex=0.1,
           layout=L$layout,vertex.size=7,vertex.color=coll[nm],vertex.label="")
      title(c("local",varname,titl))}
    
    
    if(original)
    {if (R) dev.new()
      plot(L$g2,vertex.label.color="black",vertex.label.cex=0.1,
           layout=L$layout,vertex.size=7,vertex.color=col[nm],vertex.label="")
      title(varname,titg)}
  }
  #  plot(L$layout,type="n")
  #  text(L$layout,labels=name)
  list(L,R2)
}


### getgraphnames ###
getgraphnames<-function(v,m=0,M=Inf) {
  u<-sapply(v, function(a) { 
    g<-get(a); 
    ans<-FALSE;
    if (is.igraph(g)) {
      g<-simplify(getgiant(as.undirected(g)));
      ans<-((vcount(g)<M)&(vcount(g)>m)&(assortativity.degree(g)>0))} 
    ans})
  v[u][!is.na(v[u])]
  v #several graphs
}









######################## June 16 ###########################

#### ---- read in the DataTSClean file (9/3) ---- ####
data("destrieux")
dest<-as.data.frame(destrieux)
rownames(dest)<-1:dim(dest)[1]
coords3d<-as.matrix(destrieux[,2:4])


setwd("~/Desktop/PALD")
L2data<-read.csv(file="DataTSClean2.csv",header=TRUE,row.names=1)


## Reading in each individual network

setwd("~/Desktop/PALD/BL_NET")
v<-list.files()
ML2<-NULL
L2 <-list()
for (i in 1:length(v)){
  A<-as.matrix(read_excel(v[i],col_names=FALSE))
  A<-apply(A,2,as.numeric
  rownames(A)<-1:dim(A)[1]
  colnames(A)<-1:dim(A)[1]
  ML2<-rbind(ML2,as.vector(A))
  g<-graph.adjacency(A,weighted=TRUE)
  g<-as.undirected(g) 
  g$name<-v[i]
  g$A <- A
  g$dg<-as.character(L2data$DX_bl[i])
  g$age<-as.numeric(L2data$AGE[i])
  g$gender<-as.character(L2data$PTGENDER[i])
  g$education<-as.character(L2data$PTEDUCAT[i])
  q<-as.numeric(L2data$DX_bl[i])
  g$dgn<-q
  if(q==1){g$dgn<-5}
  if(q==2){g$dgn<-1}
  if(q==3){g$dgn<-3}
  if(q==4){g$dgn<-4}
  if(q==5){g$dgn<-2}
  
 V(g)$name<-sapply(V(g)$name,function(a) gsub("V","",a))
  g$name <- i
  L2[[i]] <- g
}

i<-sample(120,1);g<-L2[[i]];cbind(t(L2data[i,c(1,4:7)]),c(g$name,g$dg,g$age,g$gender,g$education));print(g$dgn)

## look at the weighted degree

M<-t(sapply(L2,function(g) apply(g$A,1,rank)))  #this is the weight matrix

M<-t(sapply(L2,function(g) apply(g$A,1,sum)))  #this is the weight matrix

M<-t(sapply(L2,function(g) as.vector(shortest.paths(g))))

M <-t(sapply(L2,function(g) as.vector(betweenness(g))))

M <-t(sapply(L2,function(g) as.vector(get.adjacency(g,attr="weight"))))


dgn <- sapply(L2,function(g) g$dgn)
dg <- sapply(L2,function(g) g$dg)
age <- sapply(L2,function(g) g$age)
gender <- (sapply(L2,function(g) g$gender)=="Male")+1
education <- sapply(L2,function(g) g$education)

getcolour<-function(v){
  hc<-colorRampPalette(brewer.pal(8, "PiYG"))(35)
  vc<-as.numeric(v)
  col<-hc[round((vc-min(vc))/(max(vc)-min(vc))*35)+1]
  col}
Rcolors<-c(12,17,24,26,30,41,50,55,72,73,78,
           81,98,99,107,117,125,132,139,145,150,257,374,381,403,430,
           450,451,456,461,465,469,477,491,500,504,516,548, 
           554,564,578,586,596,613,619,632,638,642,657)
plot(cbind(0,1:length(Rcolors)),type="n")
text(cbind(0,1:length(Rcolors)),labels=colors()[Rcolors])
text(cbind(0,1:length(Rcolors)),labels=colors()[Rcolors],col=colors()[Rcolors])
colourv<-c("green","blue","purple","orange","red")


# only look for good and bad ones
good<- (dgn==1)|(dgn==5)
good <- good & (gender == 1)

good<-(dgn==1)
bad<-(dgn==5)
compare<-(dgn==1)|(dgn==5)

Mg<-M[good,]
Mg<-M
#rownames(Mg)<-1:120
#lay<-L$layout
Mts<-t(scale(t(Mg)))

vacross<-148*(0:73)+75+(0:73)

Macross<-Mts[,vacross]
## cluster distance
hcl_ord <- hclust(dist(Macross))$order
## an image plot of Macross
rownames(Macross) <- 1:120
myImagePlot(Macross[hcl_ord,]) 

## compare different groups in education
vector_group1 <- c(77,62,35,45,109,103,65)
myImagePlot(Macross[vector_group1,]); title(paste(vector_group1,  collapse = ", ")) 

vector_group2 <- c(77,62,35,45,109,103,65)
myImagePlot(Macross[vector_group2,]); title(paste(vector_group2, collapse = ", "))  

## compare different groups in diagnosis
dng_group1 <- c(62,35,45,109,55,103,65) #bad 1
myImagePlot(Macross[dng_group1,]); title(paste(c(dng_group1, "bad"), collapse = ", ")) 

dng_group4 <- c(62,35,45,109,55,103,65, 1, 54, 89, 2, 34) #bad 2
myImagePlot(Macross[dng_group4,]); title(paste(c(dng_group4, "bad"), collapse = ", ")) 

dng_group2 <- c(114, 102, 119, 120, 115) #good
myImagePlot(Macross[dng_group2,]); title(paste(c(dng_group2,"good"), collapse = ", ")) 

dng_group3 <- c(43, 28, 21, 31, 26) #good group2
myImagePlot(Macross[dng_group3,]); title(paste(c(dng_group3,"good"), collapse = ", ")) 

dng_group5<- which(bad)
myImagePlot(Macross[dng_group5,]); title(paste(c(dng_group5, "bad"), collapse = ", ")) 

dng_group6<- which(good)
myImagePlot(Macross[dng_group6,]); title(paste(c(dng_group6, "bad"), collapse = ", ")) 

dng_group7<- c(90, 64, 89, 1, 54, 2, 34, 32, 35, 45, 109, 103, 65, 55, 14, 29, 3, 118, 117, 19, 100, 22, 10, 15, 24)
myImagePlot(Macross[dng_group7,]); title(paste(c(dng_group5, "bad"), collapse = ", ")) 

gdiag<-L_dgn$g
A<-get.adjacency(gdiag)
v<-as.numeric(V(gdiag)$name)
A<-A[order(v),order(v)]
gdiag<-graph.adjacency(A,weighted=TRUE,mode="undirected")
D<-shortest.paths(gdiag)

## for good
a<-is.element(as.numeric(rownames(D)),which(dgn==5))

D2 <- D[a,a]

pts<-isoMDS(as.dist(D2),k=1)$points
v<-as.vector(pts)
names(v)<-rownames(pts)

## sue the order to draw our image plot
v_goodOrder <- sort(v)
v_badOrder <- c(90, 64, 89, 1, 54, 2, 34, 35, 62, 45, 109, 65, 103, 55, 14, 29, 3, 118, 117, 19, 100, 22, 10, 15, 24)
dng_groupOrdered<-as.numeric(v_badOrder) 
dng_groupOrdered<-as.numeric(names(v_goodOrder) )
myImagePlot(rbind(max(Macross),Macross[dng_groupOrdered,],min(Macross))); title(paste(c(dng_groupOrdered, "good"), collapse = ", ")) 



good<- (dgn==1)|(dgn==5)
Mg<-M[good,]
Mts<-t(scale(t(Mg)))

vacross<-148*(0:73)+75+(0:73)

Macross<-Mts[,vacross]
#Macross<-Mts
## an image plot of Macross
rownames(Macross) <- which(good)

pm<-function(a,b){par(mfrow=c(a,b))};pm(1,1)
L_dgn<-pald(cbind(Macross),dist=FALSE,scale=FALSE, show.plot=FALSE,gt=colourv[dgn[good]], main="diagnosis", show.cont=FALSE,lab=TRUE)
g2<-L_dgn$g;g3<-L_dgn$g2
A<-get.adjacency(g2,attr="weight")
v<-as.numeric(V(g2)$name)
A<-A[order(v),order(v)]
g2<-graph.adjacency(A,weighted=TRUE,mode="undirected")

A<-get.adjacency(g3,attr="weight")
v<-as.numeric(V(g3)$name)
A<-A[order(v),order(v)]
g3<-graph.adjacency(A,weighted=TRUE,mode="undirected")

lay<-L_dgn$layout
lay<-lay[order(v),]
plot(g3,vertex.label.color=colourv[dgn[good]],layout=lay,vertex.size=.8
     ,vertex.label.cex=1.5,vertex.label=which(good))
title("Figure 1: diagnosis")

gdiag<-L_dgn$g
A<-get.adjacency(gdiag)
v<-as.numeric(V(gdiag)$name)
A<-A[order(v),order(v)]
gdiag<-graph.adjacency(A,weighted=TRUE,mode="undirected")
D<-shortest.paths(gdiag)

## cluster by groups
cl<-cluster_louvain(g2)$membership
plot(g3,vertex.label.color=colourv[cl],layout=lay,vertex.size=.8,vertex.label.cex=1.5,vertex.label=which(good))


## randomly assign diagnosis
R <- NULL
v <- c(rep(1, 23), rep(0,25))
for(i in 1:1000){
  q <- sample(v)
  R<- rbind(R, c(mean(q[1:21]), mean(q[22:37]), mean(q[38:48])))
}

mean(R[,1]>=15/21)
mean(R[,2]<=3/16)
round(R, 3)
sort(R[,2])
## the result shows that they are all better than random (all lower than upper bound)

cl2<-V(g2)$name[cl==2]
cl1<-V(g2)$name[cl==1]
cl3<-V(g2)$name[cl==3]

plot(Macross[1,], ylim=c(0, 21.7), type="n")
Ma3<-Macross[is.element(rownames(Macross),cl3),]
apply(Ma3,1,lines,col="grey90")
lines(apply(Ma3,2,mean),lwd=2,col="blue")
title("Group 3")

plot(Macross[1,], ylim=c(0, 21.7), type="n")
Ma2<-Macross[is.element(rownames(Macross),cl2),]
apply(Ma2,1,lines,col="grey90")
lines(apply(Ma2,2,mean),lwd=2,col="blue")
title("Group 2")

plot(Macross[1,], ylim=c(0, 21.7), type="n")
Ma1<-Macross[is.element(rownames(Macross),cl1),]
apply(Ma1,1,lines,col="grey90")
lines(apply(Ma1,2,mean),lwd=2,col="blue")
title("Group 1")


plot(g2,vertex.label.color=colourv[dgn[good]],layout=lay,vertex.size=.8
     ,vertex.label.cex=1.2,vertex.label=education[which(good)])
title("Figure 4: diagnosis by edu")

plot(g2,vertex.label.color=colourv[3-gender[good]],layout=lay,vertex.size=.8,
     vertex.label.cex=1.5,vertex.label=which(good))
title("Figure 2: gender")

plot(g2,vertex.label.color=colourv[dgn[good]],layout=lay,vertex.size=.8
     ,vertex.label.cex=1.2,vertex.label=gender[which(good)])
title("Figure 5: diagnosis by gender")

plot(g2,vertex.label.color=colourv[dgn[good]],layout=lay,vertex.size=.8
     ,vertex.label.cex=1.25,vertex.label=round(age[which(good)]))
title("Figure 6: diagnosis by age")


edu_m <- median(as.numeric(education)[good])

edu<-as.numeric((as.numeric(education)>edu_m))+1
plot(g2,vertex.label.color=colourv[edu[good]],layout=lay,
     vertex.size=.8,vertex.label.cex=1.5,vertex.label=which(good))
title("Figure 3: education")



## this part will be for good and bad
good<-which(dgn==1)
bad<-which(dgn==5)
compare<-c(good,bad)
#compare<-1:(dim(M)[1])

## this part will be for males and females
good<-which(gender==1) # this is female
good<- good[order(dgn[good])]
bad<-which(gender==2) # this is male
bad<- bad[order(dgn[bad])]
compare<-c(good,bad)
     

Mg<-M[compare,]
for(i in 1:20)
{
  i<-8
  Mgt<-Mg[,((i-1)*148+1):(i*148)]
  L<-pald(cbind(Mgt),dist=FALSE,scale=TRUE, show.plot=FALSE,gt=colourv[dgn[compare]], 
          main=i, show.cont=FALSE)
}

setwd("~/Desktop")
pdf(file="imageplots_gender.pdf")
for (i in 1:148){
  M1<-Mg[,((i-1)*148+1):(i*148)]+1
  M1s<-t(scale(t(M1),scale=TRUE))
Mgt<-log(M1s-min(M1s)+1);myImagePlot(Mgt);title(i)
}
dev.off()



setwd("~/Desktop")
pdf(file="imageplots.pdf")
T <- NULL
for (i in 1:74){
  M1<-Mg[,((i-1)*148+1):(i*148)]+1
  M1s<-t(scale(t(M1),scale=TRUE))[,75:148]
  T <- cbind(T,apply(M1s,1,sum))
  Mgt<-log(M1s-min(M1s)+1);myImagePlot(Mgt);title(i)
}

for (i in 75:148){
  M1<-Mg[,((i-1)*148+1):(i*148)]+1
  M1s<-t(scale(t(M1),scale=TRUE))[,1:74]
  T <- cbind(T,apply(M1s,1,sum))
  Mgt<-log(M1s-min(M1s)+1);myImagePlot(Mgt);title(i)
}

Tt <-log(T-min(T)+1);myImagePlot(Tt);title("sum")
dev.off()

for(i in 1:48)
hist(T[i,],xlim=range(T),main=i,nclass=15)
## Tiz out those information !!!



## look for particular nodes to some other nodes
## we are getting connections to some particular nodes
## e.g. node 25 (area 103)

  # we want to look at detailed structure happening between nodes
#lay<-L$layout
L<-pald(cbind(Mg),dist=FALSE,scale=TRUE, show.plot=FALSE,gt=colourv[dgn[compare]], main="diagnosis", show.cont=FALSE)
lay<-L$layout
L<-pald(cbind(Mg),dist=FALSE,scale=TRUE, show.plot=FALSE,gt=colourv[gender[compare]], L=lay,main="gender", show.cont=FALSE)
L<-pald(cbind(Mg),dist=FALSE,scale=TRUE, show.plot=FALSE,gt=getcolour(education[compare]),L=lay, main="education")
L<-pald(cbind(Mg),dist=FALSE,scale=TRUE, show.plot=FALSE,gt=getcolour(age[compare]),L=lay, main="age")


ug<-apply(M[good,],2,mean)
ub<-apply(M[bad,],2,mean)

# only look for good and bad ones -- and get some commparisons
compare <- (rbind(ug-ub))
image(compare)




#image plots
Mo <- M[order(dgn),]

imageWithLabels(Mo[,1:100], row.labels = dg)
L<-pald(cbind(M),dist=FALSE,show.plot=TRUE,gt=gt2)
gt<-gt2[as.numeric(V(L$g)$name)]
hc<-colorRampPalette(brewer.pal(8, "PiYG"))(5)
plot(L$g2,vertex.color=hc[gt],vertex.label="",vertex.size=6)
plot(gt2,L$depths)

Q<-cbind(L$depths,L2data)
Q<-Q[order(L$depths),]

graph5<- hist(L$depths[gt2==5],nclass=10,xlim=c(0,1),ylim=c(0,10))
graph4<- hist(L$depths[gt2==4],nclass=10,xlim=c(0,1),ylim=c(0,10))
graph3<- hist(L$depths[gt2==3],nclass=10,xlim=c(0,1),ylim=c(0,10))
graph2<- hist(L$depths[gt2==2],nclass=10,xlim=c(0,1),ylim=c(0,10))
graph1<- hist(L$depths[gt2==1],nclass=10,xlim=c(0,1),ylim=c(0,10))

## change the color for easier visualization

plot(graph.lattice(5),vertex.color=1:5)
gt<-gt2[order(as.numeric(V(L$g)$name))]

plot(L$g2,vertex.color=hc[gt],vertex.label="",vertex.size=6)

table(gt)

## --------  want to know if same color is connected to each other ---------


E<- get.edgelist(L$g)
head(E)
E <- apply(E, 2, as.numeric)
head(E)
E2 <- cbind(E, E(L$g)$weight)
head(E2)

hist(E2[,3])
# look at the vertex of the gt value
gto2 <- gto[order(as.numeric(V(L$g2)$name))]
# larger the number, more connected
assortativity(L$g2, gto2)

# use the correct color
V(L$g2)$color <- hc[gto2]
plot(L$g2)
V(L$g)$color <- hc[gto2]
plot(L$g, vertex.label= " ", vertex.size = 5) 

# layout with short edges
plot(L$g2, vertex.label= " ", vertex.size = 5, layout=layout_with_kk(L$g2))


#--- we would like to explore more on the patterns 
#-- the depth, the severity of each disease, how are they related
#-- what does the outliers tell us ?

### get a graph for the first layout and then see the relationships ###
Layout_healthy <- layout_with_fr(L2[[1]])
### plot the graph for that layout
plot(L2[[2]],layout=Layout_healthy,vertex.size=0.1)
plot(L2[[3]],layout=Layout_healthy,vertex.size=0.1)
plot(L2[[4]],layout=Layout_healthy,vertex.size=0.1)
plot(L2[[5]],layout=Layout_healthy,vertex.size=0.1)

## want to loop through and give a type 

pdf(file = "graph types comparisons5.pdf")
for (i in which(gt2==5)){
  w <- E(L2[[i]])$weight
  w[w<quantile(w,.9)]<-0
  g<-L2[[i]]
  q<-delete.edges(g,which(w==0))
  plot(q, layout = Layout_healthy, vertex.size=0.1, vertex.label = "",edge.width=4*w/max(w)); title(gt2[i])
}

dev.off()

image2 <-myImagePlot

qu<-.9
As<-0
#pdf(file = "graph types comparisons_M1.pdf")
for (i in which(gt2==1)){
  w <- E(L2[[i]])$weight
  w[w<quantile(w,qu)]<-0
  g<-L2[[i]]
  q<-delete.edges(g,which(w==0))
  A<-get.adjacency(q,attr="weight")
  # store the adjacency matricesA into As
  As<-As+A
  print(i)
  
  # Ad is the difference
  V<-Ad>200000
  plot(image(A*V),main="1")
#  plot(q, layout = Layout_healthy, vertex.size=0.1, vertex.label = "",edge.width=4*w/max(w)); title(gt2[i])
}

As<-As/length(which(gt2==1))
image(As,main="1")
As10<-As
#dev.off()
As<-0
#pdf(file = "graph types comparisons_M1.pdf")
for (i in which(gt2==5)){
  w <- E(L2[[i]])$weight
  w[w<quantile(w,qu)]<-0
  g<-L2[[i]]
  q<-delete.edges(g,which(w==0))
  A<-get.adjacency(q,attr="weight")
  As<-As+A
  print(i)
  image2(as.matrix(A),main="1")
  #  plot(q, layout = Layout_healthy, vertex.size=0.1, vertex.label = "",edge.width=4*w/max(w)); title(gt2[i])
}
As<-As/length(which(gt2==5))
image2(as.matrix(As),main="5")
As50<-As
Ad<-(As10-As50)
#Ad[Ad<quantile(Ad,.9)]<-0
Ad<-as.matrix(Ad)
image2(Ad,main="diff")
#dev.off()


####### ------ compute centrality ------ ####### 



Lout <- lapply(L2[1:2], function(g) cent_compute(g))
L<-Lout
M<-L[[1]]
sdzero<-which(apply(M,2,function(v) sd(v)==0))
colnames(M)<-abbreviate(colnames(M),15)
M<-as.data.frame(M)
visualize_heatmap(M[,-sdzero],scale=TRUE)


L2out <- lapply(L2[-13], function(g) cent_compute(g))


####### ------ Pald Smoothing on brain network's contrality measures ------ #######
for (i in c(1:19)){
  M<-L2out[[i]]
  r<-apply(M,1,function(a) sum(is.na(a))==0)
  M<-M[r,]

  M<-apply(M,2,as.numeric)
  R2<-paldsmooth(M,V,global=FALSE,diagz = FALSE)[[2]]
  rownames(R2)<-1:dim(R2)[1]
  R3<-R2[order(R2[,5]),]
  rl<-as.numeric(rownames(R3))
  
  fileName <- gsub("xxx", i, "Sxxx_CentMeasure.pdf")
  pdf(gsub("xxx", i, "Sxxx.pdf"), fonts=c("serif", "Palatino"))
  R2<-paldsmooth(M,V[,rev(rl)],original=TRUE,diagz=FALSE)[[2]]
  dev.off()
}





####### ------ some useful functions ------ #######
#system.time({ plot(g1)}) # running time
#do.call(rbind, ListName)  # row bind elements in list
#pdf("L22.pdf", fonts=c("serif", "Palatino")); dev.off() # save plots to pdf





fastD2<-
  function(G, S = V(G),tr=FALSE,cut=3)
  {
    E(G)$weight<-1
    G<-simplify(getgiant(as.undirected(G)))
    nv <- vcount(G); nm <- length(S);
    require("Matrix", quietly = TRUE, warn.conflicts = FALSE);
    stopifnot(is.igraph(G), is.numeric(S) | class(S) == "igraph.vs");
    stopifnot(length(S) > 2, min(S) >= 1, max(S) <= nv);
    #  S <- sort(unique(as.integer(S)));
    S<-as.integer(S)
    # D <- (fastD2old(G))
    D<-shortest.paths(G)
    #D<-t(hitdist(G))+Ds
    rownames(D) <- NULL; colnames(D) <- NULL;
    if(nm < nv) {
      L <- graph.laplacian(G, sparse = TRUE) / degree(G, mode = "out", loops = FALSE);
      L[S,] <- 0; diag(L) <- 1;
      b <- Matrix(data = 0, nrow = nv, ncol = nm, sparse = TRUE);
      b[S,] <- D;
      D <- solve(L, b);
    }
    print("yep")
    D2 <- matrix(0, nrow = nm, ncol = nm);
    for(v in 1:nm) {
      nbrs <- as.integer(neighbors(G,S[v]));
      if(length(nbrs) <= 1) {
        D2[v,] <- D[nbrs,];
      } else {
        D2[v,] = apply(D[nbrs,], 2, mean);
      }
    }
    diag(D2) <- 0;
    rownames(D2)<-S
    colnames(D2)<-S
    #  return(pmin(D2,t(D2)));
    return(D2)
  }

lfr<-
  function(g) layout.fruchterman.reingold(g)

hitdist<-
  function(g)
  {
    g<-getgiant(as.undirected(g))
    P<-getP(g)
    MC <- new("markovchain",transitionMatrix = P,name = "MC")
    D <- t(meanFirstPassageTime(MC) )
    D
  }

#D<-shortest.paths(g)
#D<-hitdist(g)

g<-sawmill;V(g)$name<-1:vcount(g);lay=lfr(g);
u<-as.numeric(V(g)$gt)
V(g)$gt<-u
D<-as.matrix((fastD2(g)));  ## Community-relative distance
#D<-as.matrix(dist(D))

L<-pald(D,dist=TRUE,L=lay,global=FALSE,lab=TRUE,show.cont=FALSE)
D<-D+as.matrix(dist(D));L<-pald(D,dist=TRUE,L=lay,global=FALSE,lab=TRUE,show.cont=FALSE)
compare(L$clusters,V(g)$gt,"nmi");
plot(g,layout=lay,vertex.label.color=as.numeric(V(g)$gt)+2,vertex.size=.01,vertex.label=V(g)$name);
g2<-L$g2
A<-get.adjacency(g2)
t<-order(as.numeric(rownames(A)))
A<-A[t,t];g2<-graph.adjacency(A,mode="undirected")
D<-(fastD2(g2));
dev.new();L<-pald(pmax(D,t(D)),dist=TRUE,lab=TRUE,L=lay,show.cont=FALSE)
compare(L$clusters,V(g)$gt,"nmi");

