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

#source("functions.R")
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

#########################################################
ui <-

  navbarPage("FinnSURV-Assess PWN",
             theme = shinythemes::shinytheme("united"),

             tabPanel("1. Upload data and define parameter values",
                      withMathJax(),
                      tabsetPanel(
                        tabPanel(h5("1.1 Wood sampling"),
                                 fluidRow(h3("")),
                                 fluidRow(column(3,
                                                 wellPanel(h4(tags$b("The number of inspected sites")),
                                                           tags$hr(),
                                                           help_csv_N_w(),
                                                           fileInput("file_N_w","Upload a csv file",
                                                                     multiple = FALSE,
                                                                     accept = c("text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".csv")),
                                                           helpText("Scroll down to see the uploaded table")),

                                                 wellPanel(help_site_w(),
                                                           numericInput(inputId = "site_w",
                                                                        label = "The size of inspection site, km\\(^2\\)",
                                                                        value = 0.35, min = 0, max = 5, step = 0.1)),
                                                 wellPanel(
                                                   help_TSe_w(),
                                                   numericInput(inputId = "TSe_w",
                                                                label = "Test sensitivity",
                                                                value = 1, min = 0, max = 1, step = 0.01))),

                                          column(3,wellPanel(
                                            fluidRow(column(11,
                                                            help_n_w(),
                                                            h4(tags$b("The number of wood objects sampled per inspected site")),
                                                            tags$hr(),
                                                            radioButtons(inputId = "select_n_w_input_type",
                                                                         label = "",
                                                                         choices = c("Upload a csv file" = 1,
                                                                                     "Estimate as a probability distribution" = 2),
                                                                         selected = 2),
                                                            tags$hr(),
                                                            uiOutput(outputId = "upload_n_w_data"))))),

                                          column(3,conditionalPanel(condition = "input.select_n_w_input_type == '2'",
                                                                    wellPanel(h5(tags$i("The estimated probability distribution of the number of wood objects sampled per inspected site")),
                                                                              plotOutput(outputId = "n_wood", height=300), align = "center")))
                                 ),

                                 h5(textOutput("Table_legend_N_w")),
                                 tableOutput("table_N_w"),

                                 conditionalPanel(condition = "input.select_n_w_input_type == '1'",
                                                  h5(textOutput("Table_legend_n_w")),
                                                  tableOutput("table_n_w"))
                        ),

                        tabPanel(h5("1.2", tags$i("Monochamus"), "trapping"),
                                 fluidRow(h3("")),
                                 fluidRow(column(3,
                                                 wellPanel(h4(tags$b("The number of inspected sites")),
                                                           tags$hr(),
                                                           help_csv_N_M(),
                                                           fileInput("file_N_M", "Upload a csv file",
                                                                     multiple = FALSE,
                                                                     accept = c("text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".csv")),
                                                           helpText("Scroll down to see the uploaded table")),
                                                 wellPanel(help_site_M(),
                                                           numericInput(inputId = "site_M",
                                                                        label = "The size of inspection site, km\\(^2\\)",
                                                                        value = 0.25, min = 0, max = 5, step = 0.1)),
                                                 wellPanel(help_TSe_M(),
                                                           numericInput(inputId = "TSe_M",
                                                                        label = "Test sensitivity",
                                                                        value = 1, min = 0, max = 1, step = 0.01))
                                 ),

                                 column(3,wellPanel(
                                   fluidRow(column(11,
                                                   help_n_M(),
                                                   h4(tags$b("The number of", tags$i("Monochamus"), "sampled per inspected site")),
                                                   tags$hr(),
                                                   radioButtons(inputId = "select_n_M_input_type",
                                                                label = "",
                                                                choices = c("Upload a csv file" = 1, "Estimate as a probability distribution" = 2),
                                                                selected = 1),
                                                   tags$hr(),
                                                   uiOutput(outputId = "upload_n_M_data"))))),

                                 column(3,conditionalPanel(condition = "input.select_n_M_input_type == '2'",
                                                           wellPanel(h5(tags$i("The estimated probability distribution of the number of", tags$i("Monochamus"), "sampled per inspected site")),
                                                                     plotOutput(outputId = "n_Monochamus",height=300),align = "center")))
                                 ),

                                 h5(textOutput("Table_legend_N_M")),
                                 tableOutput("table_N_M"),

                                 conditionalPanel(condition = "input.select_n_M_input_type == '1'",
                                                  h5(textOutput("Table_legend_n_M")),
                                                  tableOutput("table_n_M"))
                        ),

                        tabPanel(h5("1.3 Target population"),
                                 fluidRow(h3("")),
                                 fluidRow(column(6,h4("At the level of inspection sites"), align = "center"),
                                          column(3,h4("At the level of regions"), align = "center")),

                                 fluidRow(column(3,
                                                 wellPanel(fluidRow(column (11,help_d_w(),
                                                                            h4(tags$b("The density of wood objects suitable for sampling")),
                                                                            tags$hr(),
                                                                            radioButtons(inputId = "select_d_w_input_type",
                                                                                         label = "",
                                                                                         choices = c("Use the estimate from Hannunen and Tuomola (2020)" = 1, "Estimate as a probability distribution" = 2),
                                                                                         selected = 1),
                                                                            tags$hr(),
                                                                            uiOutput(outputId = "upload_d_w_data"))))),

                                          column(3,
                                                 wellPanel(fluidRow(column(11,help_d_M(),
                                                                           h4(tags$b("The density of adult", tags$i("Monochamus"), "beetles")),
                                                                           tags$hr(),
                                                                           radioButtons(inputId = "select_d_M_input_type",
                                                                                        label = "",
                                                                                        choices = c("Use the estimate from Hannunen and Tuomola (2020)" = 1, "Estimate as a probability distribution" = 2),
                                                                                        selected = 1),
                                                                           tags$hr(),
                                                                           uiOutput(outputId = "upload_d_M_data"))))),


                                          column(3,
                                                 wellPanel(help_host_area(),
                                                           h4(tags$b("The area with host plants, km", tags$sup("2"))),
                                                           hr(),
                                                           help_csv_host_area(),
                                                           fileInput("file_host_area", "Upload a csv file",
                                                                     multiple = FALSE,
                                                                     accept = c("text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".csv")),
                                                           helpText("Scroll down to see the uploaded table")))
                                 ),
                                 fluidRow(column(3,
                                                 wellPanel(
                                                   h5(tags$i("The estimated probability distribution of the density of wood objects suitable for sampling")),
                                                   plotOutput("d_wood", height=300), align = "center")),
                                          column(3,
                                                 wellPanel(
                                                   h5(tags$i("The estimated probability distribution of the density of Monochamus adults")),
                                                   plotOutput("d_Monochamus", height=300), align = "center"))),

                                 fluidRow(column(12,
                                                 h5(textOutput("Table_legend_host_area")),
                                                 tableOutput("table_host_area")))
                        ),

                        tabPanel(h5("1.4 Entry sites"),
                                 fluidRow(h3("")),
                                 fluidRow(column(3,
                                                 wellPanel(help_entry_sites(),
                                                           h4(tags$b("The area of entry sites, km",tags$sup("2"))),
                                                           tags$hr(),
                                                           help_csv_entry_sites(),
                                                           fileInput("file_entry_sites", "Upload a csv file",
                                                                     multiple = FALSE,
                                                                     accept = c("text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".csv")),
                                                           helpText("Scroll down to see the uploaded table")))),
                                 fluidRow(column(12,
                                                 h5(textOutput("Table_legend_entry_sites")),
                                                 tableOutput("table_entry_sites")))
                        ))),

             tabPanel("2. Define the aim of the survey",

                      fluidRow(
                        column(3,
                               wellPanel(help_survey_type(),
                                         h4(tags$b("The aim of the surveys is")),
                                         tags$hr(),
                                         radioButtons(inputId = "select_survey_type",
                                                      label = "",
                                                      choices = c("to provide evidence to justify import requirements related to PWN and
                                                                  to facilitate export to countries with corresponding requirements" = 1,
                                                                  "to detect possible PWN invasions early enough to enable successful eradication" = 2),
                                                      selected = ""))),
                        column(3,
                               uiOutput(outputId = "localDP")),
                        column(3,
                               uiOutput(outputId = "globalDP")))),

             tabPanel("3. View results",

                      fluidRow(column(2,
                                      wellPanel(help_PriorPfree(),
                                                numericInput(inputId = "Prior_Pfree",
                                                             label = "Initial prior probability of freedom",
                                                             value = 0.5, min = 0, max = 1, step = 0.05)),
                                      wellPanel(help_Pinv(),
                                                h5(tags$b("Mean time between invasions, years")),
                                                numericInput(inputId = "Finv_min",
                                                             label = "Min",
                                                             value = 2, min = 2, max = 100, step = 1),
                                                numericInput(inputId = "Finv_max",
                                                             label = "Max",
                                                             value = 50, min = 2, max = 10000, step = 10)),
                                      wellPanel(help_n_i(),
                                                numericInput(inputId = "n_i",
                                                             label = "Number of iterations",
                                                             value = 10, min = 100, max = 10000, step = 100))),

                               column(1,
                                      actionButton(inputId = "run", label = "RUN"),
                                      tags$hr(),
                                      actionButton(inputId = "cancel", label = "CANCEL")),

                               column(8,
                                      tabsetPanel(
                                        tabPanel("Sensitivity - Country",
                                                 p(),
                                                 plotOutput(outputId = "SSe_c"),
                                                 tags$hr(),
                                                 helpText("The dots denote the medians and the bars the 95% confidence intervals of the assessment results")),
                                        tabPanel("Probability of freedom - Country",
                                                 p(),
                                                 plotOutput(outputId = "Pfree_c"),
                                                 tags$hr(),
                                                 helpText("The colored area shows the 95% confidence intervals of the assessment results")),
                                        tabPanel("Sensitivity - Regions",
                                                 p(),
                                                 plotOutput(outputId = "SSe_r", height = "1000px"),
                                                 tags$hr(),
                                                 helpText("The dots denote the medians and the bars the 95% confidence intervals of the assessment results")),
                                        tabPanel("Probability of freedom - Regions",
                                                 p(),
                                                 plotOutput(outputId = "Pfree_r", height = "1000px"),
                                                 tags$hr(),
                                                 helpText("The colored area shows the 95% confidence intervals of the assessment results"))
                                      )))),

             tabPanel("4. Download results",

                      fluidRow(column(3,
                                      wellPanel(h4(tags$b("The sensitivity of the annual surveys")),
                                                tags$hr(),
                                                selectInput("SSe_table_to_dl","2.5, 50 and 97.5% fractiles as separate csv files",
                                                            choices = c("2.5%",
                                                                        "50%",
                                                                        "97.5%")),
                                                downloadButton("download_SSe_tables", "Download"),
                                                tags$hr(),
                                                h5(tags$b("All the fractiles as a rds file")),
                                                downloadButton("download_SSe_fractiles", "Download"),
                                                tags$hr(),
                                                h5(tags$b("All iterations as a rds file")),
                                                downloadButton("download_SSe_iterations", "Download")
                                      )),
                               column(3,
                                      wellPanel(h4(tags$b("The probability of freedom after the last survey")),
                                                tags$hr(),
                                                selectInput("Pfree_table_to_dl","2.5, 50 and 97.5% fractiles as separate csv files",
                                                            choices = c("2.5%",
                                                                        "50%",
                                                                        "97.5%")),
                                                downloadButton("download_Pfree_tables", "Download"),
                                                tags$hr(),
                                                h5(tags$b("All the fractiles as a rds file")),
                                                downloadButton("download_Pfree_fractiles", "Download"),
                                                tags$hr(),
                                                h5(tags$b("All iterations as a rds file")),
                                                downloadButton("download_Pfree_iterations", "Download")
                                      )),
                               column(3,
                                      wellPanel(h4(tags$b("Download figures")),
                                                tags$hr(),
                                                help_dl_figs(),
                                                selectInput("fig_to_dl", "Choose a figure",
                                                            choices = c("Sensitivity - Country",
                                                                        "Probability of freedom - Country",
                                                                        "Sensitivity - Regions",
                                                                        "Probability of freedom - Regions")),
                                                downloadButton("download_fig", "Download"))
                               ))),

             tabPanel("About FinnSURV-Assess PWN",
                      fluidRow(column(5,
                                      tags$br(),
                                      wellPanel(
                                        p(tags$b("FinnSURV-Assess PWN is a tool for assessing the confidence in freedom from pine wood nematode
                                                 (PWN",",", tags$i("Bursaphelenchus xylophilus"),") gained in official quarantine pest surveys
                                                 in areas where PWN is not expected to cause symptoms.")),
                                        p(),
                                        p("FinnSURV-Assess PWN can be used to assess both the sensitivity of annual surveys and the probability
                                          of freedom achieved in multiannual surveys. It assumes that a) the surveys are composed of inspections
                                          that cover a fixed sized area and b) in the inspections one or more wood or Monochamus samples are collected.
                                          If all the needed data is provided separately for all regions of a country, the assessment will be done
                                          separately for each region and for the whole country."),
                                        p(),
                                        p("FinnSURV-Assess PWN was developed in the Risk Assessment Unit of the Finnish Food Authority
                                          by Salla Hannunen and Juha Tuomola in 2020. The methodology used in it is described in detail in
                                          Hannunen and Tuomola (2020) and references therein, especially",
                                          tags$a("Cannon (2002)", href = "https://doi.org/10.1016/S0167-5877(01)00262-8", target = "_blank"),",",
                                          tags$a("Martin et al. (2007)", href = "https://doi.org/10.1016/j.prevetmed.2006.09.008", target = "_blank"),",",
                                          tags$a("Efsa (2012)", href = "https://doi.org/10.2903/sp.efsa.2012.EN-366", target = "_blank"),"and",
                                          tags$a("Efsa (2018)", href = "https://doi.org/10.2903/sp.efsa.2018.EN-1399", target = "_blank"),
                                          ".")
                                        )
                                        ),
                               column(1),

                               column(5,
                                      tabsetPanel(
                                        tabPanel("Glossary",
                                                 tags$br(),
                                                 p("Apparent prevalence = The proportion of samples testing positive"),
                                                 p(),
                                                 p("Design prevalence = Roughly, design prevalence determines the minimum prevalence that the
                                                   survey is aimed to detect. If the pest prevalence is equal to or greater than the design prevalence,
                                                   at least one infested individual will be detected in the survey, with the probability equal to the
                                                   sensitivity of the survey."),
                                                 p(),
                                                 p("Early detection survey = A survey that aims to detect possible PWN invasions early enough to
                                                   enable successful eradication"),
                                                 p(),
                                                 p("Entry site = A site where the probability of PWN introduction is elevated, i.e. harbors,
                                                   industrial areas and landfills"),
                                                 p(),
                                                 p("Import-export survey = A survey that aims to provide evidence to justify import requirements
                                                   related to PWN and to facilitate export to countries with corresponding requirements"),
                                                 p(),
                                                 p("Initial prior probability of freedom = The probability that the prevalence of the pest is
                                                   below the design prevalence before the first survey"),
                                                 p(),
                                                 p("Probability of freedom = The probability that the prevalence of the pest is below the design
                                                   prevalence if the pest is not detected in the surveys"),
                                                 p(),
                                                 p("Sensitivity = Roughly, sensitivity determines the probability with which a survey is expected
                                                   to succeed in its aim. If the pest prevalence is equal to or greater than the design prevalence,
                                                   at least one infested individual will be detected in the survey, with the probability equal to
                                                   the sensitivity of the survey."),
                                                 p(),
                                                 p("Target population at the level of inspection site = Wood objects suitable for sampling / Monochamus adults"),
                                                 p(),
                                                 p("Target population at the level of regions = The area with PWN host plants in the
                                                   area for which the results of the survey will be generalized"),
                                                 p(),
                                                 p("Test sensitivity = The probability that the pest is detected in the laboratory analysis, given that it
                                                   was present in the wood object(s) / Monochamus beetle(s) included in the sample")
                                                 ),
                                        tabPanel("References",
                                                 tags$br(),
                                                 p(tags$a("Cannon RM (2002) Demonstrating disease freedom - Combining confidence levels.
                                                          Preventive Veterinary Medicine 52: 227-249.",
                                                          href = "https://doi.org/10.1016/S0167-5877(01)00262-8", target="_blank")),
                                                 p("Corine..."),
                                                 p(tags$a("European Food Safety Authority (2012)
                                                          A framework to substantiate absence of disease: the risk based estimate
                                                          of system sensitivity tool (RiBESS) using data collated according to the EFSA Standard Sample Description -
                                                          An example on", tags$i("Echinococcus multilocularis"),". Supporting Publications 2012 9(12):EN-366: 1-44.",
                                                          href = "https://doi.org/10.2903/sp.efsa.2012.EN-366", target="_blank")),
                                                 p(tags$a("European Food Safety Authority, Ciubotaru RM, Cortinas Abrahantes J, Oyedele J,
                                                          Parnell S, Schrader G, Zancanaro G, Vos S (2018) Work-plan and methodology for EFSA
                                                          to develop plant pest survey guidelines for EU Member States. EFSA supporting publication 2018 15(3):EN-1399: 1-36",
                                                          href = "https://doi.org/10.2903/sp.efsa.2018.EN-1399", target="_blank")),
                                                 p(tags$a("European Union (2012) Commission Implementing Decision of 26 September 2012 on
                                                          emergency measures to prevent the spread within the Union of",
                                                          tags$i("Bursaphelenchus xylophilus"),"(Steiner et Buhrer) Nickle et al. (the pine wood nematode)
                                                          (notified under document C(2012) 6543) (2012/535/EU).
                                                          Official Journal of the European Union L 266 2.10.2012: 42-52.",
                                                          href = "https://eur-lex.europa.eu/legal-content/EN/TXT/?uri=CELEX:02012D0535-20170310", target="_blank")),
                                                 p("Hannunen S and Tuomola J (2020)
                                                   Assessing the probability of freedom
                                                   from pine wood nematode based on 19 years of surveys. NeoBiota..."),
                                                 p(tags$a("Martin PAJ, Cameron AR, Greiner M (2007)
                                                          Demonstrating freedom from disease using multiple complex data sources
                                                          1: A new methodology based on scenario trees.
                                                          Preventive Veterinary Medicine 79: 71-97.",
                                              href = "https://doi.org/10.1016/j.prevetmed.2006.09.008", target="_blank"))
                            ),
                            tabPanel("Copyright",
                                     tags$br(),
                                     p("Please refer to FinnSURV-Assess PWN as:"),
                                     p("Hannunen S and Tuomola J 2020. FinnSURV-Assess PWN - A tool for assessing
                                       the confidence of freedom from the pine wood nematode gained in official surveys.
                                       Finnish Food Authority, Helsinki, Finland. Available at..."),
                                     tags$br(),
                                     p("FinnSURV-Assess PWN is published under", tags$a(href="https://creativecommons.org/licenses/by-nc-sa/4.0/","the Creative Commons Attribution-NonCommercial-ShareAlike 4.0
                                                                                        International", target="_blank"), "license.")
                            )
                          )
                          )))

)

#####################################################

server <- function(input,output,session){

  observe_helpers()

  # Stop session when the window is closed
  session$onSessionEnded(stopApp)

  ##########
  # INTERACTIVE PARTS OF THE UI (RENDER UI)

  # Inspection site level design prevalences of the different survey types
  output$localDP <- renderUI({
    req(input$select_survey_type)
    switch(input$select_survey_type,
           "1" = wellPanel(help_localDP_ie(),
                           h4(tags$b("Inspection site level design prevalence")),
                           tags$hr(),
                           numericInput(inputId = "DP_w",
                                        label = "Wood",
                                        value = 0.12, min = 0, max = 0.12, step = 0.01),
                           numericInput(inputId = "DP_M",
                                        label = tags$i("Monochamus"),
                                        value = 0.09, min = 0, max = 0.09, step = 0.01)),

           "2" = wellPanel(help_localDP_ed(),
                           h4(tags$b("Inspection site level design prevalence")),
                           tags$hr(),
                           numericInput(inputId = "DP_w",
                                        label = "Wood",
                                        value = 0.06, min = 0, max = 0.12, step = 0.01),
                           numericInput(inputId = "DP_M",
                                        label = tags$i("Monochamus"),
                                        value = 0.045, min = 0, max = 0.09, step = 0.01)))
  })

  # Region/country level design prevalences of the different survey types
  output$globalDP <- renderUI({
    req(input$select_survey_type)
    switch(input$select_survey_type,
           "1" = wellPanel(help_globalDP(),
                           h4(tags$b("Country level design prevalence")),
                           tags$hr(),
                           numericInput(inputId = "DPr",
                                        label = "",
                                        value = 0.01, min = 0, max = 1, step = 0.001)),
           "2" = wellPanel(help_max_inf_size(),
                           h4(tags$b("Maximum acceptable size of PWN infestation at detection, km",tags$sup("2"))),
                           tags$hr(),
                           numericInput(inputId = "max_inf_size",
                                        label = "",
                                        value = 314, min = 1, max = 5000, step = 10)))
  })

  # Number of wood objects sampled per inspected site: upload data or define as probability distribution
  output$upload_n_w_data <- renderUI({

    switch(input$select_n_w_input_type,
           "1" =  fluidRow(column(12,
                                  help_csv_n_w(),
                                  fileInput("file_n_w", "",
                                            multiple = FALSE,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv")),
                                  helpText("Scroll down to see the uploaded table"))),

           "2" =  fluidRow(
             help_pert(),
             column(6,
                    numericInput(inputId = "n_w_min",
                                 label = "Min",
                                 value = 1, min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_w_max",
                                 label = "Max",
                                 value = 5, min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_w_mode",
                                 label = "Mode",
                                 value = 2, min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_w_lambda",
                                 label = "Lambda",
                                 value = 1, min = 1, max =20, step = 1)),
             column(12,helpText("See the plot on the right for the distribution")))
    )})

  # Number of Monochamus sampled per inspected site: upload data or define as probability distribution
  output$upload_n_M_data <- renderUI({

    switch(input$select_n_M_input_type,
           "1"=    fluidRow(column(12,
                                   help_csv_n_M(),
                                   fileInput("file_n_M", "",
                                             multiple = FALSE,
                                             accept = c("text/csv",
                                                        "text/comma-separated-values,text/plain",
                                                        ".csv")),
                                   helpText("Scroll down to see the uploaded table"))),

           "2" =  fluidRow(
             help_pert(),
             column(6,
                    numericInput(inputId = "n_M_min",
                                 label = "Min",
                                 value = "", min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_M_max",
                                 label = "Max",
                                 value = "", min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_M_mode",
                                 label = "Mode",
                                 value = "", min = 1, max = 1000, step = 1)),
             column(6,
                    numericInput(inputId = "n_M_lambda",
                                 label = "Lambda",
                                 value = "", min = 1, max =20, step = 1)),
             column(12,helpText("See the plot on the right for the distribution")))
    )})

  # Density of wood objects suitable for sampling: upload data or define as porbability distribution
  output$upload_d_w_data <- renderUI({

    switch(input$select_d_w_input_type,
           "1" = fluidRow(column(12,
                                 helpText("See the plot below for the distribution"))),

           "2" = fluidRow(
             column(11,h5(tags$b("Number per km", tags$sup("2")))),
             help_pert(),
             column(6,
                    numericInput(inputId = "d_w_min",
                                 label = "Min",
                                 value = 1, min = 1, max = 10000, step = 1)),
             column(6,
                    numericInput(inputId = "d_w_max",
                                 label = "Max",
                                 value = 60, min = 1, max = 10000, step = 1)),
             column(6,
                    numericInput(inputId = "d_w_mode",
                                 label = "Mode",
                                 value = 30, min = 1, max = 10000, step = 1)),
             column(6,
                    numericInput(inputId = "d_w_lambda",
                                 label = "Lambda",
                                 value = 1, min = 1, max = 10, step = 1)),
             column(12,helpText("See the plot below for the distribution")))
    )})

  # Density of Monochamus adults: upload data or define as probability distribution
  output$upload_d_M_data <- renderUI({

    switch(input$select_d_M_input_type,
           "1"= fluidRow(column(12,
                                helpText("See the plot below for the distribution"))),

           "2" = fluidRow(column(11,h5(tags$b("Number per km", tags$sup("2")))),
                          help_pert(),
                          column(6,
                                 numericInput(inputId = "d_M_min",
                                              label = "Min",
                                              value = 1, min = 1, max = 10000, step = 1)),
                          column(6,
                                 numericInput(inputId = "d_M_max",
                                              label = "Max",
                                              value = 60, min = 1, max = 10000, step = 1)),
                          column(6,
                                 numericInput(inputId = "d_M_mode",
                                              label = "Mode",
                                              value = 30, min = 1, max = 10000, step = 1)),
                          column(6,
                                 numericInput(inputId = "d_M_lambda",
                                              label = "Lambda",
                                              value = 1, min = 1, max = 10, step = 1)),
                          column(12,helpText("See the plot below for the distribution")))
    )})

  ##########
  # UPLOADING DATA AND PRINTING IT AS TABLES

  # The number of inspected sites in the wood sampling component of the survey
  df_N_w <- reactive({
    req(input$file_N_w)
    A <- read.csv(input$file_N_w$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})

  output$Table_legend_N_w <- eventReactive(input$file_N_w,{
    "The number of inspected sites"})

  output$table_N_w <- renderTable({
    req(input$file_N_w)
    df_N_w()})

  # The number of wood objects sampled per inspected site
  df_n_w <- reactive({
    req(input$file_n_w)
    A <- read.csv(input$file_n_w$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})

  output$Table_legend_n_w <- eventReactive(input$file_n_w,{
    "The mean number of wood objects sampled per inspected site"})

  output$table_n_w <- renderTable({
    req(input$file_n_w)
    df_n_w()})

  # The number of inspected sites in the Monochamsu trapping component of the survey
  df2 <- reactive({
    req(input$file_N_M)
    A <- read.csv(input$file_N_M$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})

  output$Table_legend_N_M <- eventReactive(input$file_N_M,{
    "The number of inspected sites"})

  output$table_N_M <- renderTable({
    req(input$file_N_M)
    df2()})

  # The number of Monochamus sampled per inspected site
  df2.1 <- reactive({
    req(input$file_n_M)
    A <- read.csv(input$file_n_M$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})

  output$Table_legend_n_M <- eventReactive(input$file_n_M,{
    "The mean number of beetles sampled per inspected site"})

  output$table_n_M <- renderTable({
    req(input$file_n_M)
    df2.1()})

  # The area with host plants
  df_host_area <- reactive({
    req(input$file_host_area)
    A <- read.csv(input$file_host_area$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})

  output$Table_legend_host_area <- eventReactive(input$file_host_area,{
    "The area with host plants, km2"})

  output$table_host_area <- renderTable({
    req(input$file_host_area)
    df_host_area()})

  # The area of entry sites
  df_entry_sites <- reactive({
    req(input$file_entry_sites)
    A <- read.csv(input$file_entry_sites$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "")
    A})

  output$Table_legend_entry_sites <- eventReactive(input$file_entry_sites,{
    "The area of entry sites, km2"})

  output$table_entry_sites <- renderTable({
    req(input$file_entry_sites)
    df_entry_sites()})

  ##########
  # VECTORS NEEDED IN THE CALCULATIONS AND FIGURES

  # The considered mean times between invasions
  Finv_FI <- eventReactive(input$run,{
    seq(input$Finv_min,input$Finv_max,1)
  })

  # The years when surveys were conducted
  Y <- eventReactive(input$run,{
    a = data.matrix(df_N_w())
    seq(a[1,1],a[nrow(a),1],1)
  })

  # The number of years when surveys were conducted
  n_y <- eventReactive(input$run,{
    length(Y())
  })

  # The names of the regions included in the survey
  region_names <- eventReactive(input$run,{
    A <- read.csv(input$file_N_w$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = "'")
    colnames(A[2:ncol(A)])
  })

  # The number of regions included in the survey + 1 (to include the whole country)
  n_r <- eventReactive(input$run,{
    ncol(df_N_w())
  })

  ##########
  # UPLOADED DATA FOR THE CALCULATIONS

  # The number of inspected sites in the wood sampling component of the survey
  N_wood <- reactive({
    a = data.matrix(df_N_w())
    b = a[,2:ncol(a)]
    array(b, dim=c(nrow(b),ncol(b),input$n_i))
  })

  # The number of inspected sites in the Monochamus trapping component of the survey
  N_Monochamus <- reactive({
    a = data.matrix(df2())
    b = a[,2:ncol(a)]
    array(b, dim=c(nrow(b),ncol(b),input$n_i))
  })

  # The number of wood objects sampled per inspected site
  n_wood <- reactive({
    if(input$select_n_w_input_type == 1){
      a = data.matrix(df_n_w())
      b = a[,2:ncol(a)]
      A <- array(b, dim=c(nrow(b),ncol(b),input$n_i))
    }else{
      req(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max)
      A <- round(rpert(input$n_i,input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda))
    }
    A
  })

  # The number of Monochamus adults sampled per inspected site
  n_Monochamus <- reactive({
    if(input$select_n_M_input_type == 1){
      a = data.matrix(df2.1())
      b = a[,2:ncol(a)]
      A <- array(b, dim=c(nrow(b),ncol(b),input$n_i))
    }else{
      req(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max)
      A<- round(rpert(input$n_i,input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda))
    }
    A
  })

  # The density of wood objects suitable for sampling
  d_wood <- reactive({

    if(input$select_d_w_input_type == 1){

      # The number of Monochamus-suitable dead wood objects per km2
      obj <- round(rpert(input$n_i,166,288,398,1))

      # The proportion of Monochamus-suitable dead wood objects that is suitable for sampling
      psam <- runif(input$n_i,0.05,0.95)

      # The density of wood objects suitable for sampling, number/km2
      round(obj*psam)

    }else{
      req(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max)
      round(rpert(input$n_i,input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda))
    }
  })

  # The number of wood objects suitable for sampling per inspection site
  p_wood <- reactive({
    b <- array(d_wood(),dim=c(n_y(),n_r()-1,input$n_i))
    c = round(b*input$site_w)
    c[c < ceiling(1/input$DP_w)] <- ceiling(1/input$DP_w)
    c
  })

  # The density of Monochamus adults
  d_Monochamus <- reactive({

    if(input$select_d_M_input_type == 1){

      # The number of dead wood objects occupied by Monochamus per km2
      obju <- round(rpert(input$n_i,13.28,28.8,47.76,1))

      # The number of Monochamus eggs laid per Monochamus-suitable dead-wood object
      fobj <- rpert(input$n_i,6,31,88,1)

      # The proportion of Monochamus surviving from egg to egg-laying adults
      surv <- rpert(input$n_i,0.1,0.25,0.4,1)

      # The density of Monochamus adults, number/km2
      obju*fobj*surv

    }else{
      req(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max)
      round(rpert(input$n_i,input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda))
    }
  })

  # The number of Monochamus adults per inspection site
  p_Monochamus <- reactive({
    b <- array(d_Monochamus(),dim=c(n_y(),n_r()-1,input$n_i))
    c = round(b*input$site_M)
    c[c < ceiling(1/input$DP_M)] <- ceiling(1/input$DP_M)
    c
  })

  # The area with host plants
  host_area <- reactive({
    a = data.matrix(df_host_area())
    b = a[,2:ncol(a)]
    array(b, dim=c(nrow(b),ncol(b),input$n_i))
  })

  # The area of entry sites
  entry_sites <- reactive({
    a = data.matrix(df_entry_sites())
    b = a[,2:ncol(a)]
    matrix(b, nrow(b), ncol(b))
  })

  # The relative probability of invasion to the regions
  RP <- reactive({
    a = entry_sites()/rowSums(entry_sites())
    array(a, dim=c(nrow(a),ncol(a),input$n_i))
  })

  # The effective probability of infestation for the import-export survey
  # and the regional-level design prevalence for the early detection survey
  DPr_adj <- reactive({

    # The proportion of the target population in the regions
    PropPop = host_area()/rowSums(host_area()[,,1])

    # Weighted probability of invasion in the regions
    WR = RP()/rowSums(PropPop[,,1]*RP()[,,1])

    # Effective probability of infestation for the import-export survey
    if(input$select_survey_type == 1){
      A = input$DPr*WR
      # Regional-level design prevalence for the early detection survey
    }else{
      A = input$max_inf_size/host_area()
    }
    A
  })

  ##########
  # FIGURES FOR THE TABS "UPLOAD DATA"

  # The number of wood objects sampled per inspected site
  output$n_wood <- renderPlot({

    validate(
      need(all(input$n_w_min<=input$n_w_mode,
               input$n_w_min<=input$n_w_max,
               input$n_w_mode<=input$n_w_max,
               input$n_w_min >= 0,
               input$n_w_mode >= 0,
               input$n_w_max >= 0,
               input$n_w_lambda >= 0),
           "Make sure that min, mode, max and lambda are logical"),
      need(any(input$n_w_min!=input$n_w_mode,
               input$n_w_min!=input$n_w_max,
               input$n_w_mode!=input$n_w_max),
           "Cannot present a probability distribution if min, mode and max are the same, but analysis can be run if min, max and lambda are all defined")
    )

    x <- seq(input$n_w_min,input$n_w_max,0.1)
    dist = dpert(x,input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda)
    plot(x,dist,type = 'l',
         xlab = "Number per inspection", ylab = "Probability")
  })

  # the number of Monochamus sampled per inspected site
  output$n_Monochamus <- renderPlot({

    validate(
      need(all(input$n_M_min<=input$n_M_mode,
               input$n_M_min<=input$n_M_max,
               input$n_M_mode<=input$n_M_max,
               input$n_M_min >= 0,
               input$n_M_mode >= 0,
               input$n_M_max >= 0,
               input$n_M_lambda >= 0),
           "Make sure that min, mode, max and lambda are logical"),
      need(any(input$n_M_min!=input$n_M_mode,
               input$n_M_min!=input$n_M_max,
               input$n_M_mode!=input$n_M_max),
           "Cannot present a probability distribution if min, mode and max are the same, but analysis can be run if min, max and lambda are all defined")
    )

    x <- seq(input$n_M_min,input$n_M_max,0.1)
    dist = dpert(x,input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda)
    plot(x,dist,type = 'l',
         xlab = "Number per inspection", ylab = "Probability")
  })

  # The density of wood objects suitable for sampling
  output$d_wood <- renderPlot({
    if(input$select_d_w_input_type == 2){

      validate(
        need(all(input$d_w_min<=input$d_w_mode,
                 input$d_w_min<=input$d_w_max,
                 input$d_w_mode<=input$d_w_max,
                 input$d_w_min >= 0,
                 input$d_w_mode >= 0,
                 input$d_w_max >= 0,
                 input$d_w_lambda >= 0),
             "Make sure that min, mode, max and lambda are logical"),
        need(any(input$d_w_min!=input$d_w_mode,
                 input$d_w_min!=input$d_w_max,
                 input$d_w_mode!=input$d_w_max),
             "Cannot present a probability distribution if min, mode and max are the same, but analysis can be run if min, max and lambda are all defined")
      )

      x <- seq(input$d_w_min,input$d_w_max,0.1)
      dist = dpert(x,input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)
      plot(x,dist,type = 'l',
           xlab = expression(paste("Number per" ~ km^{2})),
           ylab = "Probability")}
    else{
      dist = round(rpert(600000,166,288,398,1))*runif(600000,0.05,0.95)
      res <- hist(dist, freq = FALSE, breaks = 80, xlim=c(0,400))
      a = res$breaks+res$breaks[2]/2
      x = a[1:length(a)-1]
      plot(x,res$density,'l',
           xlab = expression(paste("Number per" ~ km^{2})),
           ylab = "Probability", main = "")
    }
  })

  # The density of Monochamus adults
  output$d_Monochamus <- renderPlot({
    if(input$select_d_M_input_type == 2){

      validate(
        need(all(input$d_M_min<=input$d_M_mode,
                 input$d_M_min<=input$d_M_max,
                 input$d_M_mode<=input$d_M_max,
                 input$d_M_min >= 0,
                 input$d_M_mode >= 0,
                 input$d_M_max >= 0,
                 input$d_M_lambda >= 0),
             "Make sure that min, mode, max and lambda are logical"),
        need(any(input$d_M_min!=input$d_M_mode,
                 input$d_M_min!=input$d_M_max,
                 input$d_M_mode!=input$d_M_max),
             "Cannot present a probability distribution if min, mode and max are the same, but analysis can be run if min, max and lambda are all defined")
      )

      x <- seq(input$d_M_min,input$d_M_max,0.1)
      dist = dpert(x,input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)
      plot(x,dist,type = 'l',
           xlab = expression(paste("Number per" ~ km^{2})),
           ylab = "Probability")}
    else{
      dist = round(rpert(100000,13.28,28.8,47.76,1))*rpert(100000,6,31,88,1)*rpert(100000,0.1,0.25,0.4,1)
      res <- hist(dist, freq = FALSE, breaks = 50, xlim=c(0,1600))
      a = res$breaks+res$breaks[2]/2
      x = a[1:length(a)-1]
      plot(x,res$density,'l',
           xlab = expression(paste("Number per" ~ km^{2})),
           ylab = "Probability", main = "")
    }
  })

  ##########
  # CALCULATING THE SENSITIVITY OF THE ANNUAL SURVEYS

  # All iterations
  SSe_iterations <- eventReactive(input$run, {

    Sensitivity(N_wood(), n_wood(), p_wood(), input$TSe_w,
                N_Monochamus(), n_Monochamus(), p_Monochamus(), input$TSe_M,
                host_area(), RP(), DPr_adj(),
                input$select_survey_type, input$DP_w, input$DP_M, input$DPr, input$max_inf_size,
                n_r(), n_y(), input$n_i)
  })

  # 2.5, 50,97.5% fractiles
  SSe <- eventReactive(input$run, {
    Sensitivity_fractiles(SSe_iterations(),n_r(),n_y())
  })

  #####################
  # CALCULATING THE PROBABILITY OF FREEDOM AFTER THE LAST SURVEY (iterations as an Async operation)

  nclicks <- reactiveVal(0)
  Pfree_iterations_val <- reactiveVal()
  inter <- AsyncInterruptor$new()

  observeEvent(input$run,{

    req(N_wood(), n_wood(), p_wood(), input$TSe_w,
        N_Monochamus(),n_Monochamus(),p_Monochamus(),input$TSe_M,
        host_area(), entry_sites(),
        input$select_survey_type, input$DP_w, input$DP_M, any(input$DPr!="", input$max_inf_size!=""),
        input$Prior_Pfree, Finv_FI(), input$n_i)

    # Start progress monitor
    if(nclicks() == 0){
      progress <- AsyncProgress$new(session, min=1, max=100, message = "Computing")
      progress$set(message = "Computing", value = 1)
    }

    # If "Run" is clicked several times, don't do anything as analysis is already running
    if(nclicks() != 0){
      showNotification("Already running analysis")
      return(NULL)
    }
    # Increment clicks and prevent concurrent analyses
    nclicks(nclicks() + 1)

    Pfree_iterations_val(NULL)

    # Create variables that can be imported to the async operation (since reactive values and input$xxs cannot)
    N_wood <- N_wood()
    n_wood <- n_wood()
    p_wood <- p_wood()
    TSe_w = input$TSe_w
    N_Monochamus <- N_Monochamus()
    n_Monochamus <- n_Monochamus()
    p_Monochamus <- p_Monochamus()
    TSe_M = input$TSe_M
    host_area <- host_area()
    RP <- RP()
    DPr_adj <- DPr_adj()
    n_r <- n_r()
    n_y <- n_y()
    Pinv_FI <- 1/Finv_FI()
    select_survey_type = input$select_survey_type
    DP_w = input$DP_w
    DP_M = input$DP_M
    DPr = input$DPr
    max_inf_size = input$max_inf_size
    Prior_Pfree = input$Prior_Pfree
    n_i = input$n_i

    # Compute the iterations as an async operation
    Pfree_iterations <- future({

      Probability_of_freedom_iterations(progressMonitor1 = function(h) inter$execInterrupts(),
                                        progressMonitor2 = function(h) progress$set(),
                                        N_wood, n_wood, p_wood,TSe_w,
                                        N_Monochamus, n_Monochamus, p_Monochamus,TSe_M,
                                        host_area, RP, DPr_adj,
                                        select_survey_type, DP_w, DP_M, DPr, max_inf_size,
                                        Pinv_FI, Prior_Pfree,
                                        n_r, n_y, n_i)

    }) %...>% Pfree_iterations_val()

    # After the promise has been evaluated set nclicks to 0 to allow for anlother run
    # and close the progress monitor
    Pfree_iterations <- finally(Pfree_iterations,
                                function(){
                                  nclicks(0)
                                  progress$close()
                                })

    # Return something other than the promise so shiny remains responsive
    NULL
  })

  # If cancell is clicked, interrupt computing and show notification
  observeEvent(input$cancel,{
    inter$interrupt("Error: Stop Future")
    showNotification("Analysis cancelled")
  })

  # Get the 2.5, 50 and 97.5% fractiles
  Pfree <- reactive({
    req(Pfree_iterations_val())
    Pinv_FI <- 1/Finv_FI()
    Probability_of_freedom_fractiles(Pfree_iterations_val(),Pinv_FI,n_r(),n_y())
  })

  #############
  # FIGURES FOR THE TAB "VIEW RESULTS"

  # SSe - Country
  output$SSe_c <- renderPlot({

    # Validate and if not valid show an error message
    # (need(all(input$file_n_w, any(input$1, input$2,input$3)),"Error message") did not work, so I had to use if...else...)
    if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }

    req(Pfree_iterations_val())
    Plot_SSe_country(Y(),SSe(),n_r())
  })

  # Pfree - Country
  output$Pfree_c <- renderPlot({

    if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        # TAALLA
        #need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }

    req(Pfree_iterations_val())
    Plot_Pfree_country(Finv_FI(),Pfree(),n_r())
  })

  # Pfree - Regions
  output$Pfree_r <- renderPlot({

    if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }

    req(Pfree_iterations_val())
    Plot_Pfree_regions(Finv_FI(),Pfree(),region_names())
  })

  # SSe - Regions
  output$SSe_r <- renderPlot({

    if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 1 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(input$file_n_w,"Provide information on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 1){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(input$file_n_M,"Provide information on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }else if(input$select_n_w_input_type == 2 & input$select_n_M_input_type == 2){
      validate(
        need(input$file_N_w, "Upload data on the number of inspected sites in the wood sampling component"),
        need(all(input$n_w_min,input$n_w_mode,input$n_w_max,input$n_w_lambda),"Provide information on the number of wood objects sampled per inspection site"),
        need(all(input$n_w_min<=input$n_w_mode,input$n_w_min<=input$n_w_max,input$n_w_mode<=input$n_w_max),"Check the information given on the number of wood objects sampled per inspection site"),
        need(input$site_w, "Define the size of inspection site for wood sampling"),
        need(input$TSe_w, "Define test sensitivity for wood sampling"),
        need(input$file_N_M, "Upload data on the number of inspected sites in the Monochamus trapping component"),
        need(all(input$n_M_min,input$n_M_mode,input$n_M_max,input$n_M_lambda), "Provide information on the number of Monochamus sampled per inspection site"),
        need(all(input$n_M_min<=input$n_M_mode,input$n_M_min<=input$n_M_max,input$n_M_mode<=input$n_M_max),"Check the information given on the number of Monochamus sampled per inspection site"),
        need(input$site_M, "Define the size of inspection site for Monochamus trapping"),
        need(input$TSe_M,"Define test sensitivity for Monochamus trapping"),
        need(any(input$select_d_w_input_type==1,all(input$d_w_min,input$d_w_mode,input$d_w_max,input$d_w_lambda)),"Porvide information on the density of wood objects suitable for sampling"),
        need(any(input$select_d_M_input_type==1,all(input$d_M_min,input$d_M_mode,input$d_M_max,input$d_M_lambda)),"Porvide information on the density of Monochamus adults"),
        need(all(input$d_w_min<=input$d_w_mode,input$d_w_min<=input$d_w_max,input$d_w_mode<=input$d_w_max),"Check the information given on the density of wood objects suitable for sampling"),
        need(all(input$d_M_min<=input$d_M_mode,input$d_M_min<=input$d_M_max,input$d_M_mode<=input$d_M_max),"Check the information given on the density of Monochamus adults"),
        need(input$file_host_area,"Upload data on the area with host plants"),
        need(input$file_entry_sites,"Upload data on the area of entry sites"),
        need(input$select_survey_type, "Select survey type"),
        need(all(input$DP_w, input$DP_M, any(input$DPr, input$max_inf_size)),"Define design prevalences"),
        need(input$Prior_Pfree, "Determine the prior porbability of freedom"),
        need(all(input$Finv_min,input$Finv_max), "Determine the range of mean time between invasions to be considered"),
        need(input$n_i,"Determine the number of iterations")
      )
    }

    req(Pfree_iterations_val())
    Plot_SSe_regions(Y(),SSe(),region_names())
  })

  ##########
  # DOWNLOAD RESULTS

  # SSe tables
  SSetableInput <- reactive({

    a <- as.data.frame(SSe()[,,1])
    b <- as.data.frame(SSe()[,,2])
    c <- as.data.frame(SSe()[,,3])

    colnames(a)[seq(1,ncol(a)-1,1)] <- c(region_names())
    colnames(a)[ncol(a)] <- "Country"
    colnames(b) <- colnames(a)
    colnames(c) <- colnames(a)

    switch(input$SSe_table_to_dl,
           "2.5%" = a,
           "50%" = b,
           "97.5%" = c
    )
  })

  output$download_SSe_tables <- downloadHandler(
    filename = function(){
      paste(input$SSe_table_to_dl, 'csv', sep = ".")},
    content = function(file){
      write.csv(SSetableInput(), file, row.names = c(Y()))
    })

  output$download_SSe_fractiles <- downloadHandler(
    filename = function(){
      paste("SSe_fractiles.rds")},
    content = function(file){
      saveRDS(SSe(), file=file)
    })

  output$download_SSe_iterations <- downloadHandler(
    filename = function(){
      paste("SSe_iterations.rds")},
    content = function(file){
      saveRDS(SSe_iterations(), file=file)
    })

  # Pfree tables
  PfreetableInput <- reactive({

    a <- as.data.frame(Pfree()[,,1])
    b <- as.data.frame(Pfree()[,,2])
    c <- as.data.frame(Pfree()[,,3])

    colnames(a)[seq(1,ncol(a)-1,1)] <- c(region_names())
    colnames(a)[ncol(a)] <- "Country"
    colnames(b) <- colnames(a)
    colnames(c) <- colnames(a)

    switch(input$Pfree_table_to_dl,
           "2.5%" = a,
           "50%" = b,
           "97.5%" = c)
  })

  output$download_Pfree_tables <- downloadHandler(
    filename = function(){
      paste(input$Pfree_table_to_dl,"csv", sep = ".")},
    content = function(file){
      write.csv(PfreetableInput(), file, row.names = c(Finv_FI()))
    })

  output$download_Pfree_fractiles <- downloadHandler(
    filename = function(){
      paste("Pfree_fractiles.rds")},
    content = function(file){
      saveRDS(Pfree(), file=file)
    })

  output$download_Pfree_iterations <- downloadHandler(
    filename = function(){
      paste("Pfree_iterations.rds")},
    content = function(file){
      saveRDS(Pfree_iterations_val(), file=file)
    })

  # Download figures
  FigInput <- reactive({
    switch(input$fig_to_dl,
           "Sensitivity - Country" = Plot_SSe_country(Y(),SSe(),n_r()),
           "Probability of freedom - Country" = Plot_Pfree_country(Finv_FI(),Pfree(),n_r()),
           "Sensitivity - Regions" = Plot_SSe_regions(Y(),SSe(),region_names()),
           "Probability of freedom - Regions" = Plot_Pfree_regions(Finv_FI(),Pfree(),region_names())
    )})

  output$download_fig <- downloadHandler(
    filename = function(){
      paste(input$fig_to_dl,"png",sep = ".")
    },
    content = function(file){
      ggsave(FigInput(), file = file)
    })
}
############################################################
shinyApp(ui=ui,server=server)
