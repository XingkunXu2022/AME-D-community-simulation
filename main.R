if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
library(devtools)
if (!requireNamespace("amplicon", quietly=TRUE))
  install_github("microbiota/amplicon")
suppressWarnings(suppressMessages(library(amplicon)))

library(vegan)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggsci)

source('simulate_dynamics.R')
source('get_evenness.R')
source('data_processing.R')

###Specify initial community
#The abilities of the species in the simulated community to use the CH4-derived carbons (A and B)
#1:capable; 0:unable.
ability.data.frame <- data.frame(A = c(1,1,0,0,1,1),
                                 B = c(0,0,1,1,1,1))
#The average fitness of the species
fitness.data.frame <- data.frame(A = c(0.9,1,0,0,0.9,1),
                                 B = c(0,0,0.9,1,0.9,1))
#Number of the species (6)
num.species <- nrow(ability.data.frame)
#Community size (6000)
J <- num.species*1000

#Number of the generations to record
num.generation <-20

###Alpha diversity (Figure 7b)
#The data to be recorded
Pr.dif <- c()
evenness <- c()
dep <- c()
#freq.dep
freq.dep <-c(0,-2,-4,-10)

for (i in freq.dep){
  for (j in 1:1000){
    #The probability of the secreted CH4-derived carbon being A
    Pr.A <- runif(1)
    #The probability of the secreted CH4-derived carbon being B
    Pr.B <- 1-Pr.A
    
    Pr.dif_temp <- abs(Pr.A-Pr.B)
    Pr.dif <- c(Pr.dif,Pr.dif_temp)
    #Run simulation
    freq.s.list <- simulate_dynamics(num.sims = 1,  #Number of simulation
                                     num.species = num.species,  #Number of species
                                     num.generation = num.generation,  #Number of the generations to record
                                     J = J, #Community size
                                     ability.data.frame = ability.data.frame,
                                     fitness.data.frame = fitness.data.frame,
                                     freq.dep = i,
                                     Pr.c = c(Pr.A,Pr.B),
                                     rate.c = 1  #The secretion rate of the CH4-derived carbon
                                     )
    #Calculate pielou evenness index
    evenness_temp <- get_evenness_index(freq.s.list,num.generation)
    #Record evenness and dep
    evenness <- c(evenness,evenness_temp)
    dep <- c(dep,i)
    
  }
}

evenness_data <- data.frame(Pr.dif = Pr.dif,
                            evenness = evenness,
                            freq.dep = dep)
#Visualization (Figure 7b)
evenness_data$freq.dep <- as.factor(evenness_data$freq.dep)
(figure_alpha <- ggplot(data = evenness_data,aes(x=Pr.dif, y = evenness,color = freq.dep))+
  #geom_point()+
  geom_smooth()+
  labs(x=expression(paste('|',P[A],'-',P[B],'|')),
       y=paste('Pielou','\'s', 'evenness'))+
  scale_x_continuous(breaks = c(0,1.0))+
  ylim(c(0,1))+
  scale_color_lancet()+
  geom_vline(aes(xintercept=0.25),linetype='dashed')+
  geom_vline(aes(xintercept=0.5),linetype='dashed')+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
  )

ggsave('Figure 7b.png',figure_alpha,dpi=300,width = 16, height = 12,units='cm')

###Beta diversity: the relationship between the Bray-Curtis dissimilarity and |PA - PB | (Figure 7c)
#freq.dep=-4

#A matrix for output.
#Record the frequency of the species in each simulation at the 20th generation.
freq.mat <- matrix(ncol=0,nrow=num.species)

Pr.As <- c()
for (i in 1:1000){
  #The probability of the secreted CH4-derived carbon being A
  Pr.A <- runif(1)
  #The probability of the secreted CH4-derived carbon being B
  Pr.B <- 1-Pr.A
  #Record Pr.A
  Pr.As <- c(Pr.As,Pr.A)
  #Run simulation
  freq.s.list <- simulate_dynamics(num.sims = 1,  #Number of simulation
                                   num.species = num.species,  #Number of species
                                   num.generation = num.generation,  #Number of the generations to record
                                   J = J, #Community
                                   ability.data.frame = ability.data.frame,
                                   fitness.data.frame = fitness.data.frame,
                                   freq.dep = -4,
                                   Pr.c = c(Pr.A,Pr.B),
                                   rate.c = 1 #The secretion rate of the CH4-derived carbon
                                   )
  #Record the frequency of the species
  freq.mat_temp <- extract_freq.s.list(freq.s.list = freq.s.list, num.generation = num.generation)
  freq.mat <- cbind(freq.mat,freq.mat_temp)
  
}

#Calucate the Bray-Curtis dissimilarity
bray_dis <- vegdist(t(freq.mat), method = 'bray')
bray <- as.matrix(bray_dis)
#Calucate | PA1 - PA2 |
Pr.A_mat <- get_dif_matrix(values = Pr.As)
#Integrate the data of the Bray-Curtis dissimilarity and Pr.A_mat
BC_diff.Pr.A_data <- match_matrix_values(matrix1 = bray, matrix2 = Pr.A_mat)
#Visualization (Figure 7c)
(BC_dif_Pr.A <- ggplot(data = BC_diff.Pr.A_data , aes(x = matrix2, y = matrix1))+
  geom_point(color='#868686FF',alpha=1,shape=1)+
  labs(x= expression(paste('Difference in P'[A])),
       y='Bray-Curtis dissimilarity')+
  xlim(c(0,1))+
  ylim(c(0,1))+
  geom_smooth(method = 'lm', formula = y~x, se = TRUE, show.legend = FALSE,color='red') + 
  stat_poly_eq(aes(label = paste(..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y~x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top')+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
  )

ggsave('Figure 7c.png',BC_dif_Pr.A,dpi=300,width = 16, height = 12,units='cm')

#Set the secretion rate of the CH4-derived carbon
rate.c <- c(1,1.05,1.1,1.2)
#The probability of the secreted CH4-derived carbon being A
Pr.A =0.25
#The probability of the secreted CH4-derived carbon being B
Pr.B = 1-Pr.A
#Number of simulation
num.sims = 100
# A data frame for recording the frequency of the Species 1
freq.species.1 <- data.frame(freq = rep(NA,num.sims*num.generation*length(rate.c)),
                             generation = rep(NA,num.sims*num.generation*length(rate.c)),
                             rate.c = rep(NA,num.sims*num.generation*length(rate.c)))

for (i in 1:length(rate.c)){
  #Run simulation
  freq.s.list <- simulate_dynamics(num.sims = num.sims,  #Number of simulation
                                   num.species = num.species,  #Number of species
                                   num.generation = num.generation,  #Number of the generations to record
                                   J = J, #Community size
                                   ability.data.frame = ability.data.frame,
                                   fitness.data.frame = fitness.data.frame,
                                   freq.dep = -4,
                                   Pr.c = c(Pr.A,Pr.B),
                                   rate.c = rate.c[i] #The secretion rate of the CH4-derived carbon
                                   )
  #Record the frequency of the Species 1
  freq.species.1.temp <- extract.freq.species (freq.s.list, species=1)
  freq.species.1.temp$rate.c <- rep(rate.c[i],num.sims*num.generation)
  
  freq.species.1[((num.sims*num.generation)*(i-1)+1):(num.sims*num.generation*i),] <-  freq.species.1.temp
  
}

#Visualization (Figure 7d)
freq.species.1$rate.c <- as.factor(freq.species.1$rate.c)
(plot_species1 <- ggplot(data=freq.species.1,aes(x=generation,y=freq*100,color=rate.c))+
        geom_smooth()+
        ylim(0,16.7)+
        scale_x_continuous(breaks = c(1:10)*2)+
        labs(x='Generation',
             y='Frequency (%)')+
        scale_color_lancet()+
        theme_bw()+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank())
    )
ggsave('Figure 7d.png',plot_species1,dpi=300,width = 16, height = 12,units='cm')



