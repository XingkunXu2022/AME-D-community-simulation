#Reset community
reset_COM <- function(J,num.species){
  
  COM <- vector(length = J)
  for (num in 1:num.species){
    COM[((num-1)*(J/num.species)+1):(num*(J/num.species))] <- num
  }
  return(COM)
  
}

#Calculate the frequency of each species
get_freq_species <- function(COM,num.species){
  
  freq.s <- c()
  for (num in 1:num.species){
    freq.s <- c(freq.s,sum(COM == num)/length(COM))
  }
  return(freq.s)
  
}

#Output the ID of the species capable of using the secreted CH4-derived carbon
get_species_list <- function(ability.data.frame,temp.c){
  
  temp.s <- c()
  type.c <- unique(temp.c)
  for (i in type.c){
    temp.s <- c(temp.s,as.numeric(rownames(ability.data.frame)[which(ability.data.frame[,i] == 1)])) 
  }
  temp.s <- sort(unique(temp.s))
  return (temp.s)
  
}

#Calculate  probability that the individual chosen to reproduce is of one species or the other
get_Pr.s <- function(temp.c,temp.s,freq.s,freq.dep){
  
  Pr.s <- c()
  temp.freq <- freq.s[temp.s]
  temp.expectation <- c()
  
  for (num in 1:length(temp.s)){
    fit.s <- exp(freq.dep*(temp.freq[num]-0.5)+log(fitness.data.frame[temp.s[num],temp.c]))
    temp.expectation <- c(temp.expectation,fit.s*temp.freq[num])
  }
  
  for (num in 1:length(temp.s)){
    Pr <- temp.expectation[num]/sum(temp.expectation)
    Pr.s <- c(Pr.s,Pr)
  }
  
  return(Pr.s)
  
}

#The aerobic methanotrphs secrete CH4-derived carbon
secrete_compounds <- function(compound_type,rate.c,Pr.c){
  
  temp.c <- c()
  #Secrete CH4-derived carbon (A or B)
  temp.c <- c(temp.c,sample(compound_type,floor(rate.c),prob = Pr.c,replace = T))
  
  if ((rate.c - floor(rate.c)) != 0){
    if ((rate.c - floor(rate.c)) > runif(1)){
      temp.c <- c(temp.c,sample(compound_type,1,prob = Pr.c))
    }
  }
  
  return(temp.c)
  
}


#Record the frequency of the species
record_freq.s <- function(num.species,freq.s.mat,generation,COM){
  for (num in 1:num.species){
    freq.s.mat[num,generation] <- sum(COM == num)/length(COM)
  }
  return(freq.s.mat)
}


#Run simulation
simulate_dynamics <- function(num.sims,  #Number of simulation
                              num.species,  #Number of species
                              num.generation,   #Number of the generations to record
                              J, #Community size
                              ability.data.frame,  #The abilities of the species to use the CH4-derived carbons (A and B)
                              fitness.data.frame,  #The average fitness of the species
                              freq.dep = 0,  #Slpoe
                              Pr.c = NULL,  #The probability of the secreted CH4-derived carbon being A or B
                              rate.c = 1 #The secretion rate of the CH4-derived carbon
                              ){
  
  freq.s.list <- list()
  for (i in 1:num.sims){
    #Reset the matrix used for recoring the frequency of the species
    freq.s.mat <- matrix(nrow = num.species, ncol = num.generation)
    freq.s.mat[,1] <- rep(1/num.species,num.species)
    #Reset community
    COM <- reset_COM(J,num.species)
    #Reset the number of generation
    generation <- 2
    for (j in 1:(J*(num.generation-1))){
      #Calculate the frequency of each species
      freq.s <- get_freq_species(COM,num.species)
      #The aerobic methanotrphs secrete CH4-derived carbon
      temp.c <- secrete_compounds(colnames(ability.data.frame),rate.c,Pr.c)
      
      new.s <- c()
      for (compound in  temp.c){
        #Output the ID of the species capable of using the secreted CH4-derived carbon
        temp.s <- as.numeric(rownames(ability.data.frame)[which(ability.data.frame[,compound] == 1)])
        #Calculate  probability that the individual chosen to reproduce is of one species or the other
        Pr.s <- get_Pr.s(compound,temp.s,freq.s,freq.dep)
        #The ID of species chosen to reproduce 
        new.s <- c(new.s ,sample(temp.s, 1, prob = Pr.s))
      }
      
      #Select an individual to die and replace it with the individuals of the species chosen to reproduce.
      COM <- COM[-ceiling(length(COM)*runif(1))]
      COM <- c(COM,new.s)
     
      #Record the frequency of the species
      if ( j %% J == 0){
        freq.s.mat<- record_freq.s(num.species,freq.s.mat,generation,COM)
        generation <- generation+1
      }
    }
    
    cat('simulation',i,'is completed!','\n')
    
    freq.s.list[[i]] <- freq.s.mat
    
  }
  return(freq.s.list)
}

