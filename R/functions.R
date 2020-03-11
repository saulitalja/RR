
library(shiny)
library(shinythemes) # Theme
library(magrittr) # IS THIS NEEDED?
library(shinyhelper) # Pop up helps
library(ggplot2) # Plotting
library(mc2d) # Pert- probability distribution
library(ipc) # Asyc computing
library(promises) # Asyc computing
library(future) # Asyc computing
plan(multiprocess) # Asyc computing

###################################################
# A function that calculates the sensitivities of annual surveys
Sensitivity <- function(N_wood, n_wood, p_wood, TSe_wood,
                        N_Monochamus, n_Monochamus, p_Monochamus, TSe_Monochamus,
                        Pop_data, RP, DPr_adj,
                        scenario, DP_wood, DP_Monochamus, DPr, max_inf_size, 
                        n_r, n_y, n_i){
  
  #####
  # CALCULATING INSPECTION SENSITIVITY
  # Deciding whether to use binomial or hypergeometric distribution
  # (Binomial is appropriate when the sample size is less than 10% of the total population size) 
  
  if(max(n_wood)/min(p_wood) < 0.1){
    ISe_wood  = 1 - (1 - (DP_wood*TSe_wood))^n_wood
  }else{
    ISe_wood = 1 - (1 - ((n_wood*TSe_wood)/(p_wood - 0.5*(p_wood*DP_wood*TSe_wood-1))))^(p_wood*DP_wood)
  }
  
  if(max(n_Monochamus)/min(p_Monochamus) < 0.1){
    ISe_Monochamus = 1 - (1 - (DP_Monochamus*TSe_Monochamus))^n_Monochamus
  }else{
    ISe_Monochamus = 1 - (1 - ((n_Monochamus*TSe_Monochamus)/(p_Monochamus - 0.5*(p_Monochamus*DP_Monochamus*TSe_Monochamus-1))))^(p_Monochamus*DP_Monochamus)
  }
  
  #####
  # CALCULATING THE SENSITIVITY OF ANNUAL SURVEYS IN THE REGIONS
  # (Separately for wood sampling and Monochamus trapping)
  # Deciding whether to use binomial or hypergeometric distribution
  
  if(max(N_wood)/min(Pop_data) < 0.1){
    GSe_wood  = 1 - (1 - (DPr_adj*ISe_wood))^N_wood
  }else{
    GSe_wood = 1 - (1 - ((N_wood*ISe_wood)/(Pop_data - 0.5*(Pop_data*DPr_adj*ISe_wood-1))))^(Pop_data*DPr_adj)
  }
  
  if(max(N_Monochamus)/min(Pop_data) < 0.1){
    GSe_Monochamus = 1 - (1 - (DPr_adj*ISe_Monochamus))^N_Monochamus
  }else{
    GSe_Monochamus = 1 - (1 - ((N_Monochamus*ISe_Monochamus)/(Pop_data - 0.5*(Pop_data*DPr_adj*ISe_Monochamus-1))))^(Pop_data*DPr_adj)
  }
  
  #####
  # CALCULATING THE SENSITIVITY OF ANNUAL SURVEYS IN THE REGIONS AND IN FINLAND
  # (Wood sampling and Monochamus trapping combined) 
  
  # Arrays for storing the results 
  # w = wood, m = Monochamus, wm = wood and Monochamus
  GSe_wm <- array(0,dim=c(n_y,(n_r-1),n_i))
  SSe_w_FI <- array(0,dim=c(n_y,n_i))
  SSe_m_FI <- array(0,dim=c(n_y,n_i))
  SSe_wm_FI <- array(0,dim=c(n_y,n_i))
  
  for (k in 1:n_i){
    for(t in 1:n_y){
      
      # Regions
      GSe_wm[t,,k] = 1 - (1-GSe_wood[t,,k])*(1-GSe_Monochamus[t,,k])
      
      # Finland
      # Import-export
      if(scenario == 1){
        SSe_w_FI[t,k] = 1 - prod(1-GSe_wood[t,,k])
        SSe_m_FI[t,k] = 1 - prod(1-GSe_Monochamus[t,,k])
        SSe_wm_FI[t,k] = 1 - (1-SSe_w_FI[t,k])*(1-SSe_m_FI[t,k])
        # Early detection
      }else{
        SSe_wm_FI[t,k] = sum(GSe_wm[t,,k]*RP[t,,k])
      }
    }
  }
  
  #####
  # ALL THE RESULTS IN ONE ARRAY
  # Regions in columns 1-(n_r-1), Finland in column n_r
  Sensitivity_all <- array(0,dim=c(n_y,n_r,n_i))
  Sensitivity_all[,1:(n_r-1),] = GSe_wm 
  Sensitivity_all[,n_r,] = SSe_wm_FI
  return(Sensitivity_all)}


######################################################
# A function that returns the 2.5, 50 and 97.5% fractiles of the senitivity of annual surveys
Sensitivity_fractiles <- function(Sensitivity_all,n_r,n_y){
  
  # Matrixes for storing the results
  Se_2.5 <- matrix(0,n_y,n_r)
  Se_50 <- matrix(0,n_y,n_r)
  Se_97.5 <- matrix(0,n_y,n_r)
  
  for(i in 1:n_y){
    for(j in 1:n_r){
      Se_2.5[i,j] = quantile(Sensitivity_all[i,j,],0.025)
      Se_50[i,j] = quantile(Sensitivity_all[i,j,],0.5)
      Se_97.5[i,j] = quantile(Sensitivity_all[i,j,],0.975)
    }
  }
  
  # All the results in one array
  Se_all <- array(0,dim=c(n_y,n_r,3))
  Se_all[,,1] = Se_2.5 
  Se_all[,,2] = Se_50
  Se_all[,,3] = Se_97.5
  
  return(Se_all)}

#####################################################
# A function that returns the probability of freedom after the last survey
Probability_of_freedom_iterations <- function(progressMonitor1=function(h) cat("."),
                                              progressMonitor2=function(h) cat("."),
                                              N_wood, n_wood, p_wood, TSe_wood,
                                              N_Monochamus, n_Monochamus, p_Monochamus, TSe_Monochamus,
                                              Pop_data, RP, DPr_adj,
                                              scenario, DP_wood, DP_Monochamus, DPr, max_inf_size,
                                              Pinv_FI, PrioPfree_1,
                                              n_r, n_y, n_i){
  
  ##########
  # ARRAYS FOR STORING RESULTS
  # w = wood, m = Monochamus, wm = wood and Monochamus
  
  GSe_wm <-array(0,dim=c(n_y,(n_r-1),n_i))
  SSe_w_FI <-array(0,dim=c(n_y,n_i))
  SSe_m_FI <-array(0,dim=c(n_y,n_i))
  SSe_wm_FI <-array(0,dim=c(n_y,n_i))
  Pfree_w <-array(0,dim=c(n_y,(n_r-1),n_i))
  Pfree_wm <-array(0,dim=c(n_y,(n_r-1),n_i))
  Pfree_adj <-array(0,dim=c(n_y,(n_r-1),n_i))
  Pfree_w_FI <-array(0,dim=c(n_y,n_i))
  Pfree_wm_FI <-array(0,dim=c(n_y,n_i))
  Pfree_adj_FI <-array(0,dim=c(n_y,n_i))
  Pfree_all <- array(0,dim=c(length(Pinv_FI),n_r,n_i))
  
  ##########
  # CALCULATING INSPECTION SENSITIVITY
  # Deciding whether to use binomial or hypergeometric distribution
  # (Binomial is appropriate when the sample size is less than 10% of the total population size) 
  
  if(max(n_wood)/min(p_wood) < 0.1){
    ISe_wood  = 1 - (1 - (DP_wood*TSe_wood))^n_wood
  }else{
    ISe_wood = 1 - (1 - ((n_wood*TSe_wood)/(p_wood - 0.5*(p_wood*DP_wood*TSe_wood-1))))^(p_wood*DP_wood)
  }
  
  if(max(n_Monochamus)/min(p_Monochamus) < 0.1){
    ISe_Monochamus = 1 - (1 - (DP_Monochamus*TSe_Monochamus))^n_Monochamus
  }else{
    ISe_Monochamus = 1 - (1 - ((n_Monochamus*TSe_Monochamus)/(p_Monochamus - 0.5*(p_Monochamus*DP_Monochamus*TSe_Monochamus-1))))^(p_Monochamus*DP_Monochamus)
  }
  
  #####
  # CALCULATING THE SENSITIVITY OF ANNUAL SURVEYS IN THE REGIONS
  # (Separately for wood sampling and Monochamus trapping)
  # Deciding whether to use binomial or hypergeometric distribution
  
  if(max(N_wood)/min(Pop_data) < 0.1){
    GSe_wood  = 1 - (1 - (DPr_adj*ISe_wood))^N_wood
  }else{
    GSe_wood = 1 - (1 - ((N_wood*ISe_wood)/(Pop_data - 0.5*(Pop_data*DPr_adj*ISe_wood-1))))^(Pop_data*DPr_adj)
  }
  
  if(max(N_Monochamus)/min(Pop_data) < 0.1){
    GSe_Monochamus = 1 - (1 - (DPr_adj*ISe_Monochamus))^N_Monochamus
  }else{
    GSe_Monochamus = 1 - (1 - ((N_Monochamus*ISe_Monochamus)/(Pop_data - 0.5*(Pop_data*DPr_adj*ISe_Monochamus-1))))^(Pop_data*DPr_adj)
  }
  
  #########
  # CALCULATING THE SENSITIVITY OF ANNUAL SURVEYS IN THE REGIONS AND COUNTRY
  # (Wood sampling and Monochamus trapping combined)
  for (k in 1:n_i){
    for(t in 1:n_y){
      
      # Regions
      GSe_wm[t,,k] = 1 - (1-GSe_wood[t,,k])*(1-GSe_Monochamus[t,,k])
      
      # Finland
      # Import-export
      if(scenario == 1){
        SSe_w_FI[t,k] = 1 - prod(1-GSe_wood[t,,k])
        SSe_m_FI[t,k] = 1 - prod(1-GSe_Monochamus[t,,k])
        SSe_wm_FI[t,k] = 1 - (1-SSe_w_FI[t,k])*(1-SSe_m_FI[t,k])
        # Early detection
      }else{
        SSe_wm_FI[t,k] = sum(GSe_wm[t,,k]*RP[t,,k])
      }
    }  
  }
  
  # CALCULATING THE PROBABILITY OF FREEDOM
  for (h in 1:length(Pinv_FI)){   
    
    progressMonitor1(h)
    progress$set(value = (h/length(Pinv_FI))*100)
    
    Pinv_country = Pinv_FI[h]
    
    # The probability of invasion in the regions
    Pinv_regions = RP*Pinv_country
    
    for(t in 1:n_y){   
      if(t==1){  
        
        # Probability of freedom, Finland
        Pfree_wm_FI[t,] = PrioPfree_1 / (PrioPfree_1 + ((1-PrioPfree_1)*(1-SSe_wm_FI[t,])))
        
        # Adjusted probability of freedom, Finland
        Pfree_adj_FI[t,] = 1 - ( (1-Pfree_wm_FI[t,]) + Pinv_country - ((1-Pfree_wm_FI[t,])*Pinv_country) )
        
        # Probability of freedom, Regions
        Pfree_w[t,,] = PrioPfree_1 / (PrioPfree_1 + ((1-PrioPfree_1)*(1-GSe_wood[t,,])))
        Pfree_wm[t,,] = Pfree_w[t,,] / (Pfree_w[t,,] + ((1-Pfree_w[t,,])*(1-GSe_Monochamus[t,,])))
        
        # Adjusted probability of freedom, Regions
        Pfree_adj[t,,] = 1 - ( (1-Pfree_wm[t,,]) + Pinv_regions[t,,] - ((1-Pfree_wm[t,,])*Pinv_regions[t,,]) )
        
      }else{
        
        # Probability of freedom, Finland
        Pfree_wm_FI[t,] = Pfree_adj_FI[t-1,] / (Pfree_adj_FI[t-1,] + ((1-Pfree_adj_FI[t-1,])*(1-SSe_wm_FI[t,])))
        
        # Adjusted probability of freedom, Finland
        Pfree_adj_FI[t,] = 1 - ( (1-Pfree_wm_FI[t,]) + Pinv_country - ((1-Pfree_wm_FI[t,])*Pinv_country) )
        
        # Probability of freedom, Regions
        Pfree_w[t,,] = Pfree_adj[t-1,,] / (Pfree_adj[t-1,,] + ((1-Pfree_adj[t-1,,])*(1-GSe_wood [t,,])))
        Pfree_wm[t,,] = Pfree_w[t,,] / (Pfree_w[t,,] + ((1-Pfree_w[t,,])*(1-GSe_Monochamus[t,,])))
        
        # Adjusted probability of freedom, Regions
        Pfree_adj[t,,] = 1 - ( (1-Pfree_wm[t,,]) + Pinv_regions[t,,] - ((1-Pfree_wm[t,,])*Pinv_regions[t,,]) )
        
      }
    }
    
    ##########
    # RETURN ALL THE RESULTS IN ONE ARRAY
    Pfree_all[h,1:(n_r-1),] = Pfree_wm[n_y,,]
    Pfree_all[h,n_r,] = Pfree_wm_FI[n_y,]
    
  }
  
  #inter$destroy()
  #plan(sequential)
  
  return(Pfree_all)} 

######################################################
# A function that returns tha 2.5, 50 and 97.5% fractiles of the porbability of freedom after the last survey
Probability_of_freedom_fractiles <- function(Pfree_all,Pinv_FI,n_r,n_y){
  
  # Matrix for storing the results
  PF_all <- array(0,dim=c(length(Pinv_FI),n_r,3))
  
  for(i in 1:length(Pinv_FI)){
    for(j in 1:n_r){
      PF_all[i,j,1] = quantile(Pfree_all[i,j,],0.025)
      PF_all[i,j,2] = quantile(Pfree_all[i,j,],0.5)
      PF_all[i,j,3] = quantile(Pfree_all[i,j,],0.975)
    }
  }
  
  return(PF_all)}

############################################################
# FIGURES FOR THE TAB "RESULTS"

##########
# FIGURE SSe - Country
Plot_SSe_country <- function(Y,SSe,n_r){
  plot(Y,SSe[,n_r,2],
       type = "p",
       ylim=c(0.00,1.00),
       xlab = "Years",
       ylab = "Sensitivity")
  arrows(Y, SSe[,n_r,1], Y, SSe[,n_r,3], length=0.01, angle=90, code=3)
}

##########
# FIGURE Pfree - Country
Plot_Pfree_country <- function(Finv_FI,Pfree,n_r){
  plot(Finv_FI,Pfree[,n_r,1],
       type = 'l',col = 'grey',
       ylim=c(0.00,1.00),
       xlab = "Mean time between invasions, years",
       ylab = "Probability of freedom")
  lines(Finv_FI, Pfree[,n_r,3], col='grey')
  polygon(c(Finv_FI,rev(Finv_FI)), c(Pfree[,n_r,1], rev(Pfree[,n_r,3])), col='grey', border=NA)
}

##########
# FIGURE SSe - Regoins
Plot_SSe_regions <- function(Y,SSe,regions){
  
  n_regions = length(regions)  
  par(mfrow=c(ceiling(n_regions/3),3))
  par(pin=c(2,4))
  par(mar=c(2,2,2,2))
  par(mgp=c(0.5,0.75,0)) 
  par(oma=c(2.5,3,2,1)) 
  SetLine<--9
  TitleSize<- 1 
  AxisSize<- 1.2
  Xlabels<-seq(Y[1],Y[length(Y)],2)
  Ylabels<- seq(0,1,0.2)
  
  for(i in 1:n_regions){
    plot(Y, SSe[,i,2], type = "p", pch=16, frame.plot=FALSE, ylim=c(0,1), axes=FALSE, xlab="", ylab="",main=regions[i])
    arrows(Y, SSe[,i,1], Y, SSe[,i,3], length=0.01, angle=90, code=3)
    title(line=SetLine, cex.main=TitleSize)
    axis(side=1, at=Xlabels, labels=Xlabels, cex.axis=AxisSize)
    axis(side=2, at=Ylabels, labels=Ylabels, las=1, cex.axis=AxisSize)
  }
  
  mtext('Year', side=1, outer=TRUE, line=1)
  mtext('Sensitivity', side=2, outer=TRUE, line=1.5)
}

##########
# FIGURE Pfree - Regions
Plot_Pfree_regions <- function(Finv_FI,Pfree,regions){
  
  n_regions = length(regions)  
  par(mfrow=c(ceiling(n_regions/3),3))
  par(pin=c(2,4))
  par(mar=c(2,2,2,2))
  par(mgp=c(0.5,0.75,0)) 
  par(oma=c(2.5,3,2,1)) 
  SetLine<--9
  TitleSize<-1 
  AxisSize<-1.2
  Xlabels<-seq(0,Finv_FI[length(Finv_FI)],10)
  Ylabels<-c(0.00,0.20,0.40,0.60,0.80,1.00)
  
  for(i in 1:n_regions){
    plot(Finv_FI, Pfree[,i,1], type="l", col="grey", frame.plot=FALSE, 
         ylim=c(0,1), xlim=c(Finv_FI[1],Finv_FI[length(Finv_FI)]), 
         axes=FALSE, xlab="", ylab="",main=regions[i])
    lines(Finv_FI, Pfree[,i,3], col="grey")
    polygon(c(Finv_FI,rev(Finv_FI)), c(Pfree[,i,1], rev(Pfree[,i,3])), col="grey", border=NA)
    title(line=SetLine, cex.main=TitleSize)
    axis(side=1, at=Xlabels, labels=Xlabels, cex.axis=AxisSize)
    axis(side=2, at=Ylabels, labels=Ylabels, las=1, cex.axis=AxisSize)
  }
  
  mtext('The mean time between invasions to the country, years', side=1, outer=TRUE, line=1.2)
  mtext('The probability of freedom', side=2, outer=TRUE, line=1.5)
}

###############################################
# HELP TEXTS FOR THE POP UP WINDOWS

help_dl_figs <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Choose a figure",
           content = c("Note that the same figure cannot be downloaded and/or viewed twice in a row. 
                       If you need to download and/or view the same figure more than once, you need 
                       to select another figure in between the two downloads/views."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_pert <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Estimate as a probability distribution",
           content = c("If the exact value of the parameter is not known or if it varies, 
                       the parameter value can be defined as a Pert probability distribution.",
                       "",
                       "First, the minimum, maximum and mode of the parameter should be set. 
                       Then, lambda should be defined such that the resulting probability distribution 
                       reflects the uncertainty or the variation of the parameter value. Lambda 
                       takes positive values and the greater its value the more 
                       peaked is the probability distribution.",
                       "",
                       "If the exact value of the parameter is known, minimum, maximum and 
                       mode should all be set to that known value. In this special case, 
                       lambda must also be defined, but it can be given any value, 
                       as it has no effect on the outcome."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_csv_N_w <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Upload a csv file",
           content = c("Data should be uploaded as a comma separated csv file.",
                       "",
                       "In the file, data for each region should be in a separate column 
                       and data for each year should be in a separate row. The first row 
                       should have the names of the regions and the first column 
                       should have the years considered.",
                       "",
                       "All uploaded files should have the same number of columns and rows. 
                       The order of the regions should be the same in all the files 
                       and the years should always be in ascending order.",
                       "",
                       "For the regions and/or years in which no survey was done, the number of inspected 
                       sites should be zero.",
                       "",
                       "The data for the Finnish PWN surveys is in the file 
                       1.1_N_inspected_sites_wood_FI.csv."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_csv_n_w <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Upload a csv file",
           content = c("Data should be uploaded as a comma separated csv file.",
                       "",
                       "In the file, data for each region should be in a separate column 
                       and data for each year should be in a separate row. The first row 
                       should have the names of the regions and the first column 
                       should have the years considered.",
                       "",
                       "All uploaded files should have the same number of columns and rows. 
                       The order of the regions should be the same in all the files 
                       and the years should always be in ascending order.", 
                       "",
                       "For the regions and/or years in which no survey was done, the number of samples per 
                       inspected site should be zero."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_csv_N_M <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Upload a csv file",
           content = c("Data should be uploaded as a comma separated csv file.",
                       "",
                       "In the file, data for each region should be in a separate column 
                       and data for each year should be in a separate row. The first row 
                       should have the names of the regions and the first column 
                       should have the years considered.",
                       "",
                       "All uploaded files should have the same number of columns and rows. 
                       The order of the regions should be the same in all the files 
                       and the years should always be in ascending order.",
                       "",
                       "For the regions and/or years in which no survey was done, 
                       the number of inspected sites should be zero.",
                       "",
                       "The data for the Finnish PWN surveys is in the file 
                       1.2_N_inspected_sites_Monochamus_FI.csv."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_csv_n_M <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Upload a csv file",
           content = c("Data should be uploaded as a comma separated csv file.",
                       "",
                       "In the file, data for each region should be in a separate column 
                       and data for each year should be in a separate row. The first row 
                       should have the names of the regions and the first column 
                       should have the years considered.",
                       "",
                       "All uploaded files should have the same number of columns and rows. 
                       The order of the regions should be the same in all the files and 
                       the years should always be in ascending order.", 
                       "",
                       "For the regions and/or years in which no survey was done, the number of samples per 
                       inspected site should be zero.",
                       "",
                       "The data for the Finnish PWN surveys is in the file 
                       1.2_n_samples_per_inspected_site_Monochamus_FI.csv."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_csv_host_area <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Upload a csv file",
           content = c("Data should be uploaded as a comma separated csv file.",
                       "",
                       "In the file, data for each region should be in a separate column 
                       and data for each year should be in a separate row. The first row 
                       should have the names of the regions and the first column 
                       should have the years considered.",
                       "",
                       "All uploaded files should have the same number of columns and rows. 
                       The order of the regions should be the same in all the files 
                       and the years should always be in ascending order.",
                       "",
                       "Even if the the area with host plants was the same for all or some of the years, 
                       data should be given for all the considered years separately.",
                       "",
                       "The data for the Finnish PWN surveys is in the file 
                       1.3_Area_with_host_plants_FI.csv."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_csv_entry_sites <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Upload a csv file",
           content = c("Data should be uploaded as a comma separated csv file.",
                       "",
                       "In the file, data for each region should be in a separate column 
                       and data for each year should be in a separate row. The first row 
                       should have the names of the regions and the first column 
                       should have the years considered.",
                       "",
                       "All uploaded files should have the same number of columns and rows. 
                       The order of the regions should be the same in all the files 
                       and the years should always be in ascending order.",
                       "",
                       "Even if the area of entry sites was the same for all or some of the years, 
                       data should be given for all the considered years separately.",
                       "",
                       "The data for the Finnish PWN surveys is in the file 
                       1.4_Area_of_entry_sites_FI.csv."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_PriorPfree <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Initial prior probability of freedom",
           content = c("This determines the probability of freedom before the first survey.",
                       "",
                       "If no information is available about the presence/absence of PWN before 
                       the surveys were started, the initial prior probability of freedom can be set to 0.5.",
                       "",
                       "The initial prior probability of freedom should be in line with the probability of invasion, 
                       unless reason exists to assume that the probability of invasion was different before 
                       the surveys were initiated. (I.e., if the probability of invasion is assumed to be high, 
                       assuming the initial prior probability of freedom is low is not logical, and vice versa.)",
                       "",
                       "It should be noted that if the sensitivity of the annual surveys is low, even a 
                       seemingly uninformative initial prior probability of freedom (0.5) 
                       can have an impact on the probability of freedom for several years."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_Pinv <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Mean time between invasions",
           content = c("This is the mean time 
                       between events in which PWN enters the country, 
                       is transferred to a host plant and manages to establish there.",
                       "",
                       "If information or estimates about the mean time between invasions are not available, 
                       the considered range can be as wide as needed."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_n_i <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Number of iterations",
           content = c("This determines the number of iterations for the Monte Carlo simulation used in the assessment.",
                       "",
                       "For preliminary runs, a small number (<1000) is suitable. 
                       For the final assessment, 10 000 iterations are recommended.",
                       "",
                       "When the number of iterations is high, computation can take several
                       minutes, especially if a wide range of mean time between invasions is considered."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_TSe_w <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Test sensitivity",
           content = c("Test sensitivity is the probability that the pest is detected 
                       in the laboratory analysis, given that it was present in the 
                       wood object(s) included in the sample.",
                       "",
                       "Test sensitivity is always less than one. But, if inspection site level 
                       design prevalence is determined as apparent design prevalence, 
                       test sensitivity should be set to one."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_TSe_M <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Test sensitivity",
           content = c("Test sensitivity is the probability that the pest is detected 
                       in the laboratory analysis, given that it was present in the Monochamus
                       beetle(s) included in the sample.",
                       "",
                       "Test sensitivity is always less than one. But, if inspection site level 
                       design prevalence is determined as apparent design prevalence, 
                       test sensitivity should be set to one."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_site_w <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "The size of inspection site",
           content = c("This is the area covered by one inspection as square kilometers.",
                       "",
                       "The size of inspection site should be such that the number of infested 
                       wood objects per inspection site at the inspection site level design 
                       prevalence would be at least one. "),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_site_M <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "The size of inspection site",
           content = c("This is the area covered by one inspection as square kilometers.",
                       "",
                       "The size of inspection site should be such that the number of infested 
                       Monochamus adults per inspection site at the inspection site level design 
                       prevalence would be at least one. "),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_n_w <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "The number of wood objects sampled per inspected site",
           content = c("This is the number of wood objects from which material was taken for 
                       laboratory analysis per inspected site, irrespective of whether 
                       the material from all the objects was pooled or analyzed separately."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_n_M <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "The number of Monochamus sampled per inspected site",
           content = c("This is the number of Monochamus collected for laboratory analysis 
                       per inspected site, irrespective of whether all the beetles 
                       were pooled or analyzed separately."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_entry_sites <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "The area of entry sites",
           content = c("This is the area of sites in with elevated probability of PWN introduction, 
                       i.e. harbors, industrial areas and landfills, in square kilometers."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_host_area <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "The area with host plants",
           content = c("This is the area with PWN host plants within the area for which the results 
                       of the surveys will be generalized, in square kilometers."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_d_w <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "The density of wood objects suitable for sampling",
           content = c("This is the number of wood objects suitable for 
                       sampling per square kilometer in the area for which 
                       the results of the survey will be generalized."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_d_M <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "The density of adult Monochamus beetles",
           content = c("This is the number of adult Monochamus beetles per 
                       square kilometer in the area for which the results 
                       of the survey will be generalized."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_survey_type <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "The aim of the surveys",
           content = c("The first option is referred to as import-export survey and the latter as early detection survey. 
                       The pest populations that these two survey types aim to detect differ in two respects.",
                       "",
                       "First, for import-export surveys, PWN infestation is expected to be randomly distributed throughout 
                       the country, while for early detection surveys, PWN infestation is expected to be restricted to one region.",
                       "",
                       "Second, the design prevalences of import-export surveys can be based merely on what is required 
                       to facilitate trade, while those of early detection surveys should 
                       be defined based on the size of an eradicable PWN population."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_localDP_ie <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Inspection site level design prevalence",
           content = c("Inspection site level design prevalence is determined as the proportion 
                       of PWN-infested wood objects of all wood objects suitable for sampling, 
                       and the proportion of PWN-infested Monochamus adults of all Monochamus adults.",
                       "",
                       "Design prevalence should be such that PWN could reach it, at least at some point 
                       in time, if it was established in the considered area. Additionally, it 
                       must be such that it corresponds to, at least, one whole infested wood object 
                       and Monochamus adult per inspection site.",
                       "",
                       "For import-export surveys, inspection site level design prevalences can 
                       correspond to the prevalence that established PWN populations would likely have in the considered 
                       conditions.",
                       "",
                       "For areas where PWN is not expected to cause symptoms and where Bursaphelenchus mucronatus is 
                       established, inspection site level design prevalence may be set equal to observed 
                       prevalence of B. mucronatus.",
                       "",
                       "If the prevalence of B. mucronatus has been estimated with similar sampling as that done in 
                       the PWN surveys, inspection site level design prevalence can be determined as apparent prevalence 
                       and test sensitivity (on the tab 'Upload data and define parameter values') can be set to one."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_localDP_ed <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Inspection site level design prevalence",
           content = c("Inspection site level design prevalence is determined as the proportion 
                       of PWN-infested wood objects of all wood objects suitable for sampling, 
                       and the proportion of PWN-infested Monochamus adults of all Monochamus adults.",
                       "",
                       "Design prevalence should be such that PWN could reach it, at least at some point 
                       in time, if it was established in the considered area. Additionally, it 
                       must be such that it corresponds to, at least, one whole infested wood object 
                       and Monochamus adult per inspection site.",
                       "",
                       "For the early detection survey, inspection site level design prevalences should 
                       preferably correspond to a PWN population that is still growing 
                       (i.e., has not reached its maximum prevalence).",
                       "",
                       "In areas where PWN is not expected to cause symptoms and where Bursaphelenchus mucronatus 
                       is established, inspection site level design prevalence should be lower than the prevalence of B. mucronatus.",
                       "",
                       "If the prevalence of B. mucronatus has been estimated with similar sampling as that done in 
                       the PWN surveys, inspection site level design prevalence can be determined as apparent prevalence 
                       and test sensitivity (on the tab 'Upload data and define parameter values') can be set to one."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_globalDP <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Country level design prevalence",
           content = c("For import-export surveys, country level design prevalence is determined as the proportion of 
                       PWN-infested area of the total area with PWN host plants in the area 
                       for which the results of the surveys will be generalized to.",
                       "",
                       "Country level design prevalence can be based on requirements of the trading partners, 
                       political considerations, availability of resources, and/or biological plausibility.",
                       "",
                       "Country level design prevalence must be such that PWN could reach it, at least at some 
                       point in time, if it was established in the country. Also, it must be such that it corresponds 
                       to, at least, the size of one whole infested inspection site."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}

help_max_inf_size <- function(){
  tags$i()%>%
    helper(type="inline",
           size = "m",
           title = "Maximum acceptable size of PWN infestation at detection",
           content = c("For early detection surveys, region level design prevalence is defined based on 
                       the maximum size of an eradicable PWN infestation in square kilometers.",
                       "",
                       "The maximum area from which eradication could be attempted is limited by the financial 
                       and physical resources available for the eradication measures (e.g., the availability of 
                       labor and machinery). If the maximum area from which eradication could be attempted can 
                       be estimated, it may be used as a proxy for the maximum area from which eradication could be 
                       achieved. However, it is important to keep in mind that attempting and achieving are very different 
                       things when it comes to eradicating invasive species.",
                       "",
                       "In EU countries, the requirements of the EU emergency measures for PWN (EU 2012) can be used to guide the definition of 
                       the maximum size of an eradicable PWN population. The measures allow EU member states to refrain from 
                       attempting eradication of PWN if the diameter of the infested area exceeds 20 km. Hence, it is 
                       logical to assume that in an early detection survey, PWN should be detected before the diameter of 
                       its infestation is more than 20 km.",
                       "",
                       "The area of a PWN infestation that has a 20 km diameter depends on the shape of the area and the proportion 
                       of the area that is covered by PWN host plants. If the infested area is circular, it is, at maximum, 
                       314 km2. In most cases, it is much smaller due to low coverage of PWN host plants. Again, it is important 
                       to remember that it is not at all clear if such a large infestation could be eradicated with the resources 
                       available for delimiting the infested area and conducting the eradication measures."),
           buttonLabel = "OK", 
           easyClose = TRUE,
           icon = "question-circle",
           colour = "orangered",
           fade = FALSE)
}
