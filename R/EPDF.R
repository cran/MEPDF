
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

  count.coord<-plyr::ddply(coord,lst,nrow)
  count.coord$V1<-count.coord$V1/n/prod(grid.sizes)


  hst<-function(y){

    eval<-function(x){
      y<-floor(x/grid.sizes)

      prob<-count.coord[which(apply(count.coord[,1:dim],1,function(z) all(z == y))==TRUE),]$V1

      if(pracma::isempty(prob)){return(0)}
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


pseudokernel<-function(data,mn,mx,grid.sizes,rings){

  # Basic data values, dimension etc.
  dim <- ncol(data)
  n <- nrow(data)

  # Add corners if necessary
  if (missing(mx)) {
    max.corner <- rep(max(data), dim)
  }
  if (missing(mn)) {
    min.corner <- rep(min(data), dim)
  }


  # Placement function, assigning coordinates to input values
  pl<-function(x){
    return(ceiling((x-mn)/grid.sizes))
  }

  pts<-t(apply(X = data,MARGIN = 1,FUN = pl))
  gridl<-(mx-mn)/grid.sizes
  pgrid<-array(0,dim = (mx-mn)/grid.sizes)

  # Due to the array implementation, it is necessary to transform coordinates to counting vectors
  vtoc<-function(v){
    if(length(v) == 1){return(v[1])}else{
      return(v[1] + sum((v[-1]-1)*cumprod(rev(gridl)[-1])))
    }
  }

  # Add the weight on the respective grid cell
  addval<-function(v){
    pgrid[vtoc(v)]<<-pgrid[vtoc(v)] + 1/(prod(grid.sizes)*n)
  }

  invisible(apply(X = pts,MARGIN = 1,FUN = addval))


  # Superimpose layers upon which the weight is distributed
  if(!missing(rings)){

    pgrid2<-array(0,dim = (mx-mn)/grid.sizes)

    addval2<-function(v){
      tryCatch({

        p<-pgrid[vtoc(v)]/(rings+1)


        v1<-gtools::permutations(v = c(-1,0,1),r = dim,n = 3,repeats.allowed = TRUE)
        v1<-v1[-(nrow(v1)+1)/2,]
        vr<-gtools::permutations(v = seq(-rings,rings,by = 1),r = dim,n = 2*rings+1,repeats.allowed = TRUE)
        vr<-vr[-(nrow(vr)+1)/2,]

        for(r in 1:rings){
          if(r == 1){
            vmat<-v1 + v[col(v1)]
            vmat<-vmat[apply(vmat, 1, function(row) all(row > 0 )),]
            vmat<-vmat[apply(vmat, 1, function(row) all(row < rev(gridl)+1)),]

            pos<-apply(X = vmat,MARGIN = 1,FUN = vtoc)
            pgrid2[pos[intersect(pos > 0,pos < (length(pgrid2)+1))]]<<-pgrid2[pos[intersect(pos > 0,pos < (length(pgrid2)+1))]] + p/((2*r+1)^dim - (2*r-1)^dim)
          }else{
            vmat<-vr + v[col(vr)]
            vmat<-vmat[apply(vmat, 1, function(row) all(row > 0 )),]
            vmat<-vmat[apply(vmat, 1, function(row) all(row < rev(gridl)+1)),]

            pos<-apply(X = vmat,MARGIN = 1,FUN = vtoc)
            pgrid2[pos[intersect(pos > 0,pos < (length(pgrid2)+1))]]<<-p/((2*r+1)^dim - (2*r-1)^dim)
            pgrid2[pos[intersect(pos > 0,pos < (length(pgrid2)+1))]]<<-pgrid2[pos[intersect(pos > 0,pos < (length(pgrid2)+1))]] + p/((2*r+1)^dim - (2*r-1)^dim)
          }
        }

        pgrid2[vtoc(v)] <<- p
      })
    }

    invisible(apply(X = pts,MARGIN = 1,FUN = addval2))
  }


  # Return values for the standard grid ()
  if(missing(rings)){
    pdf<-function(x){
      return(pgrid[vtoc(pl(x))])
    }


    return(list(pdf = pdf,grid = pgrid))
  }else{
    pdf1<-function(x){
      return(pgrid[vtoc(pl(x))])
    }

    pdf2<-function(x){
      return(pgrid2[vtoc(pl(x))])
    }

    return(list(pdf1 = pdf1,pdf2 = pdf2,grid1 = pgrid,grid2 = pgrid2))
  }



}




normkernel<-function(x,H){

  d<-length(x)
  if(missing(H)){H<-diag(d)}

  return(mvtnorm::dmvnorm(x,mean = rep(0,d),sigma = H))
}

epakernel<-function(x,H){

  k<-function(x){
    dim<-length(x)
    ret<-(1-x^2)

    if(all(abs(x) < 1)){
      return(prod(ret)*(3/4)^dim)
    }else{return(0)}
  }


  a<-det(H)^(-1/2)
  b<-pracma::inv(H)^(1/2)

  return(a*k(b*x))
}




ekde<-function(x,data,H,rule,kernel = normkernel){

  dim<-ncol(data)
  n<-nrow(data)


  if(missing(H)){
    if(rule == "silverman"){
      H = apply(X = data,FUN = stats::var,MARGIN = 2)*diag(dim)*(4/(dim+2))^(2/(dim+4))*n^(-2/(dim+4))
    }
    if(rule == "scott"){
      H = apply(X = data,FUN = stats::var,MARGIN = 2)*diag(dim)*n^(-2/(dim+4))
    }
    if(missing(rule)){
      H = diag(dim)
    }
  }

  return(mean(apply(X = x-data,FUN = kernel,MARGIN = 1,H = H)))

}


