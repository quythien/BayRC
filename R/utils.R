adjust.to.2pi = function(x){
  d = x/(2*pi)
  d.abs = floor(abs(d))
  if(d>=0&d<1){
    return(x)
  }else if(d<0){
    return(x+(d.abs+1)*2*pi)
  }else if(d>1){
    return(x-(d.abs)*2*pi)
  }
}
