extract_freq.s.list <- function(freq.s.list,num.generation){
  
  num.species <- nrow(freq.s.list[[1]])
  extracted_freq.s <- matrix(nrow = num.species,
                              ncol = length(freq.s.list))
  i <- 1
  for (j in freq.s.list){
    extracted_freq.s[,i] <- j[,num.generation]
    i <- i+1
  }
  return(extracted_freq.s)
  
}

extract.freq.species <- function(freq.s.list, species=1){
 
  freq <- c()
  generation <- c()
  num.generation <- ncol(freq.s.list[[1]])
  
  for (i in 1:length(freq.s.list)){
    for (j in 1:num.generation ){
      freq <- c(freq,freq.s.list[[i]][species,j])
      generation <- c(generation,j)
    }
  }
  
  extracted_mat <- data.frame(freq = freq,
                              generation = generation)
  return(extracted_mat)

}

#Calculate pielou's evenness
get_evenness_index <- function(freq.s.list,num.generation){
  
  library(vegan)
  extracted_freq.s <- extract_freq.s.list(freq.s.list,num.generation)
  shannon_index <- diversity(t(extracted_freq.s),index='shannon',base=exp(1))
  richness <- nrow(freq.s.list[[1]])
  pielou <- shannon_index/log(richness,exp(1))
  return(pielou)
  
}



