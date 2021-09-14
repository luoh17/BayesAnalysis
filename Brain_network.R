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

### ----------------------------- functions ---------------------- ###

## random graph with same number of nodes (sample from all possible edges)
randgraph<-function (n,e=n) {
  g2<-random.graph.game(n,e,type="gnm")
  g2<-induced.subgraph(g2,largest_comp(g2))
  g2
}

## gives us the biggest component of the graph
largest_comp<-function(graph) {
  cl <- clusters(graph)
  V(graph)[which.max(cl$csize) == cl$membership]
}

## get the transition matrix

getP<-function (g) {
  A<-as.matrix(get.adjacency(g))
  r<-apply(A,1,sum)
  P<-diag(1/r)%*%A
  P
}

## pald 

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

## pald-smooth (runs pald pn N matrix, smooth the overlaid variables according to cohesion)
 # averaging over the pald neighbors

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

## pm produces an a x b grid of plots
pm<-function(a,b){par(mfrow=c(a,b))};
pm(1,1)


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

# minmax scaling
minmax <- function(v){
  if(sd(v) == 0){
    ans<-rep(0.5, length(v))
  }
  else{
    ans <-(v-min(v))/(max(v)-min(v))
  }
  ans
}


# Define a function for plotting a matrix 
myImagePlot<-
  function(x, gt=1, cex.axis=1,...){
    gt <- rev(gt)
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
    axis(LEFT <-2, at=1:length(yLabels), labels=FALSE, las= HORIZONTAL<-1,
         cex.axis=cex.axis)
    Map(axis, side=2, at=1:length(yLabels), col.axis=gt, labels=yLabels, lwd=0, las=1,cex.axis=cex.axis)
    
    # Color Scale
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    
    layout(1)
  }

# getcolor does a color ramp for v (from dark purple[low] to dark green[high])
# plot(1:10, col=getcolour(1:10), pch=17, cex=2)
getcolour<-function(v){
  hc<-colorRampPalette(brewer.pal(8, "PiYG"))(35)
  vc<-as.numeric(v)
  col<-hc[round((vc-min(vc))/(max(vc)-min(vc))*35)+1]
  col}
Rcolors<-c(12,17,24,26,30,41,50,55,72,73,78,
           81,98,99,107,117,125,132,139,145,150,257,374,381,403,430,
           450,451,456,461,465,469,477,491,500,504,516,548, 
           554,564,578,586,596,613,619,632,638,642,657)
pm(1,1)
plot(cbind(0,1:length(Rcolors)),type="n")
text(cbind(0,1:length(Rcolors)),labels=colors()[Rcolors])
text(cbind(0,1:length(Rcolors)),labels=colors()[Rcolors],col=colors()[Rcolors])
# color options we can use
colourv<-c("green","blue","purple","orange","red")


### ------------------------ centrality computation -------------------- ###

# g1 can be used to see all the centralities
g1<- randgraph(10)
# some of the centrality measures are excluded 
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
  
  rownames(M) <- paste(c(rownames(M),g1$name),collapse=".")
  M <- as.matrix(M)
  M
}


#### ----------------- read in the DataTSClean file -------------------- ####
data("destrieux")
# dest is all the information about the brain regions (locations)
dest<-as.data.frame(destrieux)
rownames(dest)<-1:dim(dest)[1]
coords3d<-as.matrix(destrieux[,2:4])   #3d coordinates
texts3d(coords3d,texts=dest$name,cex=.5,col=as.numeric(as.factor(dest$lobe)))


### ----- read in data --------###
setwd("~/Desktop/PALD")
# L2data is the patient information
L2data<-read.csv(file="DataTSClean2.csv",header=TRUE,row.names=1)


## -------------- Reading in each individual network (L2) --------------- ##

setwd("~/Desktop/PALD/BL_NET")
v<-list.files()

### L2 is a 120 long list of graphs with patients' information (age, weights, ...)
L2 <-list()
ML2 <- NULL
for (i in 1:length(v)){
  A<-as.matrix(read_excel(v[i],col_names=FALSE))
  A<-apply(A,2,as.numeric)
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

## get a sense of the brain structure
par(mfrow=c(2,2))
hist(log(E(L2[[1]])$weight+1), main="individual 1 data", xlab="log(edge weight +1)")
hist(log(E(L2[[2]])$weight+1), main="individual 2 data", xlab="log(edge weight +1)")
hist(log(E(L2[[3]])$weight+1), main="individual 3 data", xlab="log(edge weight +1)")
hist(log(E(L2[[4]])$weight+1), main="individual 4 data", xlab="log(edge weight +1)")

par(mfrow=c(5,5))
for (i in 1:20){
hist(log(E(L2[[i]])$weight+1), main=i, xlab="log(edge weight +1)")
}

g<-L2[[1]]
E(g)$weight[E(g)$weight<exp(12)]<-0

g3<-g
V(g3)$label<-""


open3d()
rglplot(g3,layout=coords3d,labels=dest,
        vertex.color=as.numeric(as.factor(dest$lobe)),
        edge.width=.2*log(E(g3)$weight+1))

movie3d(spin3d(axis = c(0, 0, 1)), duration = 15,
        dir = getwd())


rglplot(g3,layout=coords3d,labels=dest,
        vertex.color=as.numeric(as.factor(dest$hemi)),edge.width=E(g3))

library(magick)



# check to assure everything is working
i<-sample(120,1);g<-L2[[i]];cbind(t(L2data[i,c(1,4:7)]),c(g$name,g$dg,g$age,g$gender,g$education));print(g$dgn)


## M is the matrix, its rows represent different graphs

M<-t(sapply(L2,function(g) apply(g$A,1,rank)))  #this is the weight matrix

M<-t(sapply(L2,function(g) apply(g$A,1,sum)))  #this is the weight matrix

M<-t(sapply(L2,function(g) as.vector(shortest.paths(g))))

M <-t(sapply(L2,function(g) as.vector(betweenness(g))))

### this is what we focused on: 120 by 21904(148*148)
M <-t(sapply(L2,function(g) as.vector(get.adjacency(g,attr="weight"))))


#### vectors of patient's characteristics
# 1: cog normal 2: SMC 3: EMCI 4: LMCI 5: AD
dgn <- sapply(L2,function(g) g$dgn)
dg <- sapply(L2,function(g) g$dg)
# range 55 - 90.3
age <- sapply(L2,function(g) g$age)
# males = 1, females = 2
gender <- (sapply(L2,function(g) g$gender)=="Male")+1
# range 11 to 20
education <- sapply(L2,function(g) as.numeric(g$education))

good<- (dgn==1)|(dgn==5) # subset of individuals with particular characteristics 
#good <- good & (gender == 1)
#good<-(dgn==1)

#bad<-(dgn==5)
#compare<-(dgn==1)|(dgn==5)

Mg<-M[good,] ## Mg is the graph data only for particular patients(good)
#Mg<-M

Mgl<-log(Mg+1)
Mts<-t(scale(t(Mgl))) #Mts is the scaled graph information with the # of sd
#Mts <- t(apply(Mgl,1,minmax)) #Mts is the scaled graph information with minmax (smalles:0, biggest:1)

# for 148 x 148 data, vacross is the positions of the same region cross hemisphere weights
vacross<-148*(0:73)+75+(0:73)
# randome pairning as compared to the original

assmat<-NULL
for (i in 1:400){
#vacross <- 148*(0:73)+75+sample(0:73)
vacross<-148*(0:73)+75+(0:73)
  vacross.s<-sample((74),74)
#  vacross.s<-vacross
# Macross is Mts restricted to only cross ties (74 columns)
#vacross<-assmat[300,-1][1:3]
  Macross<-Mts[,vacross[vacross.s]]
#Macross<-t(apply(Macross, 1, minmax))
# dendrogram for Macross data
rownames(Macross) <- which(good)
#hcl <- hclust(dist(Macross), method = "complete")
#hcl_ord<-hcl$order
#plot(hcl)
#myImagePlot(Macross[hcl_ord,],col.lab="green") 
#myImagePlot(Macross[order(dgn[good]),],gt=sort(dgn[good])) 

## L_dgn is a pald list for Macross data
#L_dgn<-pald(cbind(Macross),dist=FALSE,scale=FALSE, show.plot=FALSE,gt=colourv[dgn[as.numeric(rownames(Macross))]], main="diagnosis", show.cont=FALSE,lab=TRUE)
L_dgn<-pald(cbind(Macross),dist=FALSE,scale=FALSE, show.plot=FALSE, main="diagnosis", show.cont=FALSE,lab=TRUE)
# g1 is the cohesion graph, g2 is threshold
g1<-L_dgn$g;g2<-L_dgn$g2
A<-get.adjacency(g1,attr="weight")
v<-as.numeric(V(g1)$name)
A<-A[order(v),order(v)]
g1<-graph.adjacency(A,weighted=TRUE,mode="undirected")

A<-get.adjacency(g2,attr="weight")
v<-as.numeric(V(g2)$name)
A<-A[order(v),order(v)]
g2<-graph.adjacency(A,weighted=TRUE,mode="undirected")

lay<-L_dgn$layout
lay<-lay[order(v),]
#plot(g2,vertex.label.color=colourv[dgn[good]],layout=lay,vertex.size=.8
 #    ,vertex.label.cex=1.5,vertex.label=which(good))
ass<-assortativity(g2, dgn[good])
assmat<-rbind(assmat,c(ass,vacross.s))
title(round(ass,2))
print(c(i,ass))
}
assmat<-assmat[order(assmat[,1]),];tail(assmat)

n<-dim(assmat)[1]
vacross.s<-assmat[n,-1]
vacross.s<-1:74
Macross<-Mts[,vacross[vacross.s]]
Macross<-t(apply(Macross, 1, minmax))
# dendrogram for Macross data
rownames(Macross) <- which(good)
#hcl <- hclust(dist(Macross), method = "complete")
#hcl_ord<-hcl$order
#plot(hcl)
#myImagePlot(Macross[hcl_ord,],col.lab="green") 
myImagePlot(Macross[order(dgn[good]),],gt=sort(dgn[good])) 

## L_dgn is a pald list for Macross data
#L_dgn<-pald(cbind(Macross),dist=FALSE,scale=FALSE, show.plot=FALSE,gt=colourv[dgn[as.numeric(rownames(Macross))]], main="diagnosis", show.cont=FALSE,lab=TRUE)
L_dgn<-pald(cbind(Macross),dist=FALSE,scale=FALSE, show.plot=FALSE, main="diagnosis", show.cont=FALSE,lab=TRUE)
# g1 is the cohesion graph, g2 is threshold
g1<-L_dgn$g;g2<-L_dgn$g2
A<-get.adjacency(g1,attr="weight")
v<-as.numeric(V(g1)$name)
A<-A[order(v),order(v)]
g1<-graph.adjacency(A,weighted=TRUE,mode="undirected")

A<-get.adjacency(g2,attr="weight")
v<-as.numeric(V(g2)$name)
A<-A[order(v),order(v)]
g2<-graph.adjacency(A,weighted=TRUE,mode="undirected")

lay<-L_dgn$layout
lay<-lay[order(v),]
plot(g2,vertex.label.color=colourv[dgn[good]],layout=lay,vertex.size=.8
     ,vertex.label.cex=1.5,vertex.label=which(good))
ass<-assortativity(g2, dgn[good])
#assmat<-rbind(assmat,c(ass,vacross.s))
destrieux[vacross.s,]

####################################################

lobe<-dest$lobe


rownames(Mts) <- which(good)
runn<-function(n,z)
{
  bestval<-(-1)
  bestvec<-1:z
  assvec<-NULL
for (i in 1:n){
vacross <- 148*(0:73)+75+sample(0:73)

vacross <- 148*(0:73)+75+(0:73)%%74

vacross<-sample(vacross, z)
#vacross<-assvec[300,3:(2+z)]
# Macross is Mts restricted to only cross ties (74 columns)
Macross<-Mts[,vacross]
#Macross<-Mts
Macross<-scale(Macross)
Macross<-t(apply(Macross, 1, minmax))
colnames(Macross)<-vacross
# dendrogram for Macross data
rownames(Macross) <- which(good)
hcl <- hclust(dist(Macross), method = "complete")
hcl_ord<-hcl$order
#plot(hcl)
a<-dgn[which(good)]
#dev.new();myImagePlot(Macross[order(dgn[good]),],gt=rev(sort(a))) 
#dev.new();myImagePlot(Macross[hcl_ord,],gt=rev(a[hcl_ord])) 


## L_dgn is a pald list for Macross data
L_dgn<-pald(cbind(Macross),dist=FALSE,scale=FALSE, show.plot=FALSE,gt=colourv[dgn[good]], 
show.cont=FALSE,lab=TRUE)
# g1 is the cohesion graph, g2 is threshold
g1<-L_dgn$g;g2<-L_dgn$g2
A<-get.adjacency(g1,attr="weight")
v<-as.numeric(V(g1)$name)
A<-A[order(v),order(v)]
g1<-graph.adjacency(A,weighted=TRUE,mode="undirected")

A<-get.adjacency(g2,attr="weight")
v<-as.numeric(V(g2)$name)
A<-A[order(v),order(v)]
g2<-graph.adjacency(A,weighted=TRUE,mode="undirected")

lay<-L_dgn$layout
lay<-lay[order(v),]
#plot(g2,vertex.label.color=colourv[dgn[good]],layout=lay,vertex.size=.8
#     ,vertex.label.cex=1.5,vertex.label=which(good))
#title("Figure 1: diagnosis")

assort<-assortativity(g2, dgn[good])
if (assort>bestval)
{
  gt<-dgn[good]
  pald(cbind(Macross),dist=FALSE,scale=FALSE, show.plot=TRUE,gt=colourv[dgn[good]], 
       show.cont=FALSE,lab=TRUE)
  bestval<-assort
  title(assort)
  }

destname<-dest$name[(vacross-74)%%148]

assvec<-rbind(assvec, c(assort,mean(Macross),sort(vacross),sort((vacross-74)%%148)))
print(round(sort(assvec[,1]),2))
}
assvec<-assvec[order(assvec[,1]),]
#list(assvec=assvec, vacross=vacross,destname=destname)
#}


assvec[order(assvec[,1]),]
destname<-dest$name[assvec[dim(assvec)[1],(z+3):(2*z+2)]]

list(assvec=assvec,destname=destname)
}


# ------- image plots for different groups
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

exp(10.5)



