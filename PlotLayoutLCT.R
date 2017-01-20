plot.layout.lct = function(Result,txtsize=0.5, label=Result[[1]]){
  library(diagram)
  
  result = Result[[3]][1:(length(Result[[3]])-1)]
  N_class = length(c(unlist(result),result[[length(result)]]))
  M = matrix("0",N_class,N_class)
  
  Ch = (1:result[[1]]) + 1 #which classes in the hierarchy have to be tested
  teller = 1
  index = list()
  for(i in 1:length(result)){
    for(j in 1:length(result[[i]])){
      if(i!=1){
        Ch = (max(index[[teller-1]])+1)
        length(Ch) = result[[i]][j]
        if(any(is.na(Ch))){
          for(k in which(is.na(Ch))){Ch[k] = Ch[k-1]+1}}}
      index[[teller]] = Ch
      teller = teller+1
    }}
  for(i in 1:length(index)){
    M[index[[i]],i] = ""
  }
  plotmat(M,pos=sapply(label,length),box.type="rect", box.size=0.00,arr.length=0.000001,lcol="red",name=unlist(label), box.cex=txtsize, lwd=5)
}