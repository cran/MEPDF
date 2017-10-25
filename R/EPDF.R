

cube<-function(data,mx,mn,grid.sizes){

  if(missing(mx)){mx<-max(data)}
  if(missing(mn)){mn<-min(data)}

  dim<-ncol(data)
  n<-nrow(data)
  if(missing(grid.sizes)){grid.sizes<-rep(0.05,dim)}

  allupr<-function(x){all(x < mx)}
  alllwr<-function(x){all(x > mn)}

  dat<-data[which(lapply(as.list(data.frame(t(data))),allupr) == TRUE),]
  d<-dat[which(lapply(as.list(data.frame(t(dat))),alllwr) == TRUE),]

  cube.inc<-ceiling((mx-mn)/grid.sizes)
  cube.length<-cube.inc*grid.sizes

  fill<-rep(0,prod(cube.inc))
  grid<-as.list(fill)
  dim(grid)=c(cube.inc)

  coord<-data.frame(floor(d/grid.sizes))

  # Create List argument for ddply
  char<-"X1"
  if(dim > 1){
    for(i in 2:dim){
      char<-paste(char,",X",i,sep="")
    }
  }

  lst<-eval(parse(text=paste(".(",char,")",sep="")))

  count.coord<-ddply(coord,lst,nrow)
  count.coord$V1<-count.coord$V1/n/prod(grid.sizes)


  hst<-function(y){

    eval<-function(x){
      y<-floor(x/grid.sizes)

      prob<-count.coord[which(apply(count.coord[,1:dim],1,function(z) all(z == y))==TRUE),]$V1

      if(isempty(prob)){return(0)}
      else{
        return(count.coord[which(apply(count.coord[,1:dim],1,function(z) all(z == y))==TRUE),]$V1)
      }
    }

    return(unlist(lapply(y,eval)))
  }

  return(hst)
}


epdf<-function(data,max.corner,min.corner,main.gridsize,rescubes){

  dim<-ncol(data)
  if(missing(max.corner)){max.corner<-rep(max(data),dim)}
  if(missing(min.corner)){min.corner<-rep(min(data),dim)}

  # initial grid
  main.grid<-cube(data,max.corner,min.corner,main.gridsize)

  if(missing(rescubes))
    {hst<-main.grid}
  else{
    # superimpose areas of higher resolutions
    n<-length(rescubes)
    grid.sizes<-lapply(rescubes, `[[`, 3)
    mxs<-lapply(rescubes, `[[`, 2)
    mns<-lapply(rescubes, `[[`, 1)

    for(i in 1:n){
      assign(paste("subcubes",i,sep=""),cube(data,mn=mns[[i]],mx=mxs[[i]],grid.sizes=grid.sizes[[i]]))
    }

    # Partial EPDF selection

    hst<-function(x){

      xmax<-function(y){all(unlist(x)<y)}
      xmin<-function(y){all(unlist(x)>y)}

      upr<-which(unlist(lapply(mxs,xmax)))
      lwr<-which(unlist(lapply(mns,xmin)))

      common<-intersect(lwr,upr)

      # return cube value
      if(length(common) == 0){
        return(main.grid(x))
      }
      else{
        return(eval(parse(text=paste("subcubes",max(common),sep="")))(x))
      }
    }
  }

return(hst)
}

