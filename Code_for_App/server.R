
server <- function(input, output) {
  Route     = isolate(input$route)
  Dose      = isolate(input$doselevel)
  Tdoses    = isolate(input$numdose)
  tinterval = isolate(input$tinterval)
  end_time  = isolate(input$simu_time)

  observeEvent(input$reset,{reset("Treatment")})
  
  ## Load model
  mod <- mcode ("NPsPBPK.code", NPsPBPK.code) # compile the basic NPs model
  pars <- c(BW = 0.273)
  
  ## Set up simulation subjects
  r1 <- reactive({
    input$action |input$recalc
    Tdoses    = input$numdose
    tinterval = input$tinterval
    end_time  = input$simu_time
    BW        = 0.273 # The body weight adopted from experiment
    
    # Refine the dose regiment
    Dose      = isolate (input$doselevel*BW)
    Total_Dose     = Dose*Tdoses
    end_time  = isolate (input$simu_time)
    
    dt        = 0.2
    N = isolate(input$n_sim)
    
    Route     = isolate (input$route)
    Size      = isolate (input$size_oral|input$size_IV|input$size_IT)
    ZETA      = isolate (input$Zeta)
    HD        = isolate (input$HD)
    SA        = isolate (input$SA)
    nNPs      = isolate (10^input$nNPs)
    logSize   = log(Size)
    logZETA   = log(abs(ZETA))
    logHD     = log(HD)
    logSA     = log(SA)
    lognNPs   = log(nNPs)
    
    #///////////////////////////////////////////////////////////////////////////
    # Oral exposure route
    if (Route == "Oral") {
      KurineC_eq         = -7.76e-05+3.83e-05*SA+0.000167*logHD
      KSRESrelease_eq    = 11.8+0.916*logSA-0.44*lognNPs
      KSRESmax_eq        = -52.7-4.91*logSA+2.32*lognNPs
      KRRESrelease_eq    = 3.43-0.772*logHD-0.432*logSA
      KRRESmax_eq        = 0.358+0.00822*Size-0.0173*SA
      KpulRESrelease_eq  = 19.2+0.0683*Size-5.3*logHD
      KpulRESmax_eq      = 4.2+0.0454*Size+4.91*logSA
      KLRESrelease_eq    = -0.237+0.872*SA-1.84*logSA
      KLRESmax_eq        = -41.8+1.04*SA+13.2*logZETA
      KKRESrelease_eq    = -0.0288-0.0207*ZETA-1.4e-14*nNPs
      KKRESmax_eq        = 19.5+1.03*SA-0.777*lognNPs
      KinRESrelease_eq   = 16.2+0.0224*HD-3.81*logHD
      KinRESmax_eq       = -2.2+0.108*HD+1.21*SA
      KGIb_eq            = 0.00244-0.000544*logHD-0.000234*logSA
      KfecesC_eq         = -3.93-0.01*Size+1.74*logZETA
      KbileC_eq          = 0.325+0.00226*HD+0.0143*ZETA
      
      pars = c(
        AGIREScap  = 100,
        AinREScap  = 100,
        AKREScap   = 100,
        ALREScap   = 100,
        ApulREScap = 100,
        ARREScap   = 100,
        ASREScap   = 100,
        PAGIC      = 0.001,
        PAKC       = 0.001,
        PALC       = 0.001,
        PALuC      = 0.001,
        PARC       = 0.001,
        PASC       = 0.001,
        KGIRESmax  = 1,
        KGIRESrelease = 0.001,
        Kpulrestra = 0.000865,        
        KGIb       = ifelse(KGIb_eq<0, 0.0004, KGIb_eq),
        KbileC     = ifelse(KbileC_eq<0, 0.05, KbileC_eq),
        KfecesC    = ifelse(KfecesC_eq<0, 1.32, KfecesC_eq),
        KinRESmax  = ifelse(KinRESmax_eq<0, 20, KinRESmax_eq),
        KinRESrelease = ifelse(KinRESrelease_eq<0, 7, KinRESrelease_eq),  
        KRRESmax   = ifelse(KRRESmax_eq<0, 0.1, KRRESmax_eq),
        KRRESrelease = ifelse(KRRESrelease_eq<0, 0.2, KRRESrelease_eq), 
        KurineC    = ifelse(KurineC_eq<0, 0.001, KurineC_eq)
      )
      
      ppars = c(
        KKRESmax   = ifelse(KKRESmax_eq<0, 13.51, KKRESmax_eq),
        KKRESrelease = ifelse(KKRESrelease_eq<0, 0.099, KKRESrelease_eq), 
        KLRESmax   = ifelse(KLRESmax_eq<0, 16.58, KLRESmax_eq),
        KLRESrelease = ifelse(KLRESrelease_eq<0, 9.59, KLRESrelease_eq), 
        KpulRESmax = ifelse(KpulRESmax_eq<0, 18.336, KpulRESmax_eq),
        KpulRESrelease = ifelse(KpulRESrelease_eq<0, 6.21, KpulRESrelease_eq),
        KSRESmax   = ifelse(KSRESmax_eq<0, 4.39, KSRESmax_eq),
        KSRESrelease = ifelse(KSRESrelease_eq<0, 0.89, KSRESrelease_eq) 
      )
      
      
      
      ex <- ev(ID = 1:N, amt = Dose, ii = tinterval, addl = Tdoses-1, 
               cmt = "ALumen")  
      
    }
    
    #////////////////////////////////////////////////////////////// 
    # IV route
    if (Route == "IV") {
      ## Define the multiple linear equations for each of the variables
      PASC_eq           =  1+0.00848*SA-0.0316*lognNPs
      KurineC_eq        =  1.18e-05-6.42e-06*SA+5.15e-18*nNPs
      KSRESrelease_eq   = -0.22-0.00991*ZETA+0.00407*SA
      KSRESmax_eq       = 597+31.6*logHD-206*logZETA
      KRRESrelease_eq   = -0.106+0.00292*HD+0.0309*SA
      KRRESmax_eq       = -5.67-1.23e-14*nNPs+0.26*lognNPs
      KpulRESrelease_eq = -4.17+0.0584*ZETA+1.78*logZETA
      KpulRESmax_eq     = -2.63+7.86e-15*nNPs+1.12*logHD
      KLRESrelease_eq   = 0.07-0.000305*Size-0.0211*logSA
      KLRESmax_eq       = 317+4.69*SA-10.7*lognNPs
      KKRESrelease_eq   = 4.49-0.0272*Size-1.59*logSA
      KKRESmax_eq       = 2.88+2.93e-14*nNPs-0.11*lognNPs
      KinRESrelease_eq  = -1.55-0.000414*Size+0.52*logZETA
      KinRESmax_eq      = 3.42+0.0194*Size-0.101*lognNPs
      KGIRESrelease_eq  =-1.17-0.0522*ZETA+0.0153*SA
      KGIRESmax_eq      = 48.2-0.0748*HD-1.47*lognNPs
      KGIb_eq           = 0.00541-1.19e-20*SA+9.16e-33*nNPs
      KfecesC_eq        = -1.36+0.00216*HD+0.0536*lognNPs
      KbileC_eq         = 0.0302+0.000144*Size+0.00142*ZETA
      ASREScap_eq       = -76.1+111*SA+398*logHD
      ARREScap_eq       = 5.06+3.01*SA+3.36e-13*nNPs
      ApulREScap_eq     = -448+24.7*logHD+15.7*lognNPs
      ALREScap_eq       = 8.21e+03+1.08e+03*logSA-251*lognNPs
      AinREScap_eq      = -448+24.7*logHD+15.7*lognNPs
      AGIREScap_eq      = 7.95+2.92*SA+3.25e-13*nNPs
      
      pars = c(
        AGIREScap       = ifelse(AGIREScap_eq<0, 100, AGIREScap_eq),
        AinREScap       = ifelse(AinREScap_eq<0, 100, AinREScap_eq),
        AKREScap        = 10,
        ALREScap        = ifelse(ALREScap_eq<0, 4000, ALREScap_eq),
        ApulREScap      = ifelse(ApulREScap_eq<0, 100, ApulREScap_eq),
        ARREScap        = ifelse(ARREScap_eq<0, 100, ARREScap_eq),
        ASREScap        = ifelse(ASREScap_eq<0, 4000, ASREScap_eq),
        PAGIC           = 0.001,
        PAKC            = 0.001,
        PALC            = 0.001,
        PALuC       = 0.001,
        PARC        = 0.001,
        PASC        = ifelse(PASC_eq<0, 0.25, PASC_eq),
        KGIRESmax   = ifelse(KGIRESmax_eq<0, 0.5, KGIRESmax_eq),
        Kpulrestra  = 0.000865,        
        KGIb        = ifelse(KGIb_eq<0, 5.41e-3, KGIb_eq),
        KbileC      = ifelse(KbileC_eq<0, 8.0e-4, KbileC_eq),
        KinRESmax   = ifelse(KinRESmax_eq<0, 0.2, KinRESmax_eq),
        KRRESmax    = ifelse(KRRESmax_eq<0, 2, KRRESmax_eq), 
        KKRESmax    = ifelse(KKRESmax_eq<0, 0.46, KKRESmax_eq),
        KLRESmax    = ifelse(KLRESmax_eq<0, 113, KLRESmax_eq),
        KpulRESmax  = ifelse(KpulRESmax_eq<0, 1.94e-1, KpulRESmax_eq),
        KSRESmax    = ifelse(KSRESmax_eq<0, 2, KSRESmax_eq),  
        KurineC     = ifelse(KurineC_eq<0, 1.48e-5, KurineC_eq),
        KfecesC     = ifelse(KfecesC_eq<0, 0.373, KfecesC_eq)
      )
      
      ppars = c(
        KGIRESrelease  = ifelse(KGIRESrelease_eq<0, 0.35, KGIRESrelease_eq),
        KKRESrelease   = ifelse(KKRESrelease_eq<0, 9.9e-3, KKRESrelease_eq),
        KinRESrelease  = ifelse(KinRESrelease_eq<0, 0.03, KinRESrelease_eq),
        KLRESrelease   = ifelse(KLRESrelease_eq<0, 1.09e-2, KLRESrelease_eq),
        KRRESrelease   = ifelse(KRRESrelease_eq<0, 2.45e-2, KRRESrelease_eq),
        KpulRESrelease = ifelse(KpulRESrelease_eq<0, 2.45e-2, KpulRESrelease_eq),
        KSRESrelease   = ifelse(KSRESrelease_eq<0, 0.8, KSRESrelease_eq)
      )
      
      ex <- ev(ID = 1:N, amt = Dose, ii = tinterval, addl = Tdoses-1,  
               cmt = "AV")  
      
    }
    #////////////////////////////////////////////////////////////// 
    # IT exposure route
    if (Route == "IT") {
      KurineC_eq           = 1.27e-05-1.16e-05*SA+1.02e-17*nNPs
      KSRESrelease_eq      = 0.463-0.00981*SA-4.18e-15*nNPs
      KSRESmax_eq          = 0.013-0.00035*logHD-0.000113*lognNPs
      KpulRESrelease_eq    = 0.00859-0.000338*SA+2.2e-17*nNPs
      KpulRESmax_eq        = -2.16e+03+229*logHD+60.6*lognNPs
      KKRESrelease_eq      = 0.393-0.0167*SA+0.0174*logHD
      KKRESmax_eq          = 0.424-0.692*SA+5.3e-13*nNPs
      KinRESrelease_eq     = -0.0112-0.000997*ZETA-0.000328*SA
      KinRESmax_eq         = -55-54.6*Size+53.7*HD
      KfecesC_eq           = 0.000262-1.73e-08*Size-9.72e-07*logSA
      KbileC_eq            = 0.0953+0.000741*Size+0.00469*ZETA
      ApulREScap_eq        = -458+1.17e+03*logHD+1.46e+03*logSA
      AinREScap_eq         = -458+1.17e+03*logHD+1.46e+03*logSA
      
      pars = c(
        PALuC         = 0.001, 
        PAGIC         = 0.001, 
        PALC          = 0.001, 
        PAKC          = 0.001, 
        PASC          = 0.001, 
        PARC          = 0.001, 
        AGIREScap     = 10, 
        ALREScap      = 10, 
        ASREScap      = 10, 
        AKREScap      = 10, 
        ApulREScap    = ifelse(ApulREScap_eq<0, 6500, ApulREScap_eq), 
        ARREScap      = 1.5,
        KGIb          = 1e-4, 
        KGIRESmax     = 1, 
        KinRESmax     = 300, 
        AinREScap     = ifelse(AinREScap_eq<0, 6500, AinREScap_eq), 
        KRRESmax      = 0.01, 
        KbileC        = ifelse(KbileC_eq<0, 0.0008, KbileC_eq),
        KfecesC       = ifelse(KfecesC_eq<0, 0.00026, KfecesC_eq),
        KurineC       = ifelse(KurineC_eq<0, 0.0000426, KurineC_eq), 
        Kpulrestra    = 0.05,        
        KLRESmax      = 0.00898,        
        KSRESmax      = ifelse(KSRESmax_eq<0, 0.00869, KSRESmax_eq), 
        KKRESmax      = ifelse(KKRESmax_eq<0, 0.0449, KKRESmax_eq), 
        KpulRESmax    = ifelse(KpulRESmax_eq<0, 278, KpulRESmax_eq)
        
      )
      
      ppars <-c(
        KLRESrelease   = 0.0706,         
        KRRESrelease   = 0.5, 
        KSRESrelease   = ifelse(KSRESrelease_eq<0, 0.13, KSRESrelease_eq), 
        KKRESrelease   = ifelse(KKRESrelease_eq<0, 0.0706, KKRESrelease_eq), 
        KinRESrelease  = 0.002,        
        KGIRESrelease  = 0.1, 
        KpulRESrelease = ifelse(KpulRESrelease_eq<0, 0.00178, KpulRESrelease_eq) 
      )
      
      ev_1 <- ev(ID = 1:N, amt = 0.54*Dose, ii = tinterval, addl = Tdoses-1, cmt = "Atra")
      ev_2 <- ev(ID = 1:N, amt = 0.46*Dose, ii = tinterval, addl = Tdoses-1, cmt = "Apul")
      ex <- ev_1 + ev_2 
      
    }
    #/////////////////////////////////////////////////////////////
    # Inhalation exposure route
    if (Route == "IH") {
      Size = isolate (input$size_IH)
      if (Size == "23 nm") {
        pars = c(
          PALuC = 0.001, PAGIC = 0.001, PALC = 0.001, PAKC = 0.001, 
          PASC = 0.001, PARC = 0.001, ALREScap = 100, ASREScap = 100, 
          AKREScap = 100, KpulRESrelease = 0.1, KpulRESmax = 300, ApulREScap = 5000, KinRESrelease = 0.1, KinRESmax = 300, 
          AinREScap = 5000, KRRESrelease = 0.1, KRRESmax = 1, ARREScap = 1.5,
          AGIREScap = 100, KbileC = 3e-4,  KurineC = 5e-3,
          KfecesC = 0.132, Kpulrestra = 0.003)
        
        ppars <- c(KLRESrelease = 0.004, KLRESmax = 0.07,  
                   KSRESrelease = 0.001, KSRESmax = 1, 
                   KKRESrelease = 0.0001, KKRESmax = 0.07, 
                   KGIRESrelease = 0.5, KGIRESmax = 0.1,
                   KGIb = 1e-4)
      }
      ev_1 <- ev(ID = 1:N, amt = 0, ii = tinterval, addl = Tdoses-1, cmt = "Aua") 
      ev_2 <- ev(ID = 1:N, amt = 0.06822*Dose, ii = tinterval, addl = Tdoses-1, cmt = "Atra")
      ev_3 <- ev(ID = 1:N, amt = 0.08178*Dose, ii = tinterval, addl = Tdoses-1, cmt = "Apul")
      ex <- ev_1 + ev_2 + ev_3
    } 
    
    #////////////////////////////////////////////////////////
    ## Simulation
    sd.log.pars <- sqrt(log(1 + (ppars*0.5)^2/ ppars^2))
    m.log.pars  <- log(ppars) - 0.5*sd.log.pars^2
    
    ##////////////////////////////////////////////////////////
    a    <- vector()
    idata <- vector()
    idata <- tibble(ID=1:N)
    for( i in 1:length(names(ppars))) {
      a = 
        rlnormTrunc(
          N,
          meanlog = m.log.pars,
          sdlog   = sd.log.pars,
          min = qlnorm(0.025, meanlog = m.log.pars, sdlog = sd.log.pars),
          max = qlnorm(0.975, meanlog = m.log.pars, sdlog = sd.log.pars)
        )
      
      idata <- cbind.data.frame(idata, a)
    }
    
    colnames(idata)[2:length(names(idata))] <- c(names(ppars))    
    #////////////////////////////////////////////////////////////
    tsamp = tgrid(0, tinterval*(Tdoses-1) + end_time*24, dt)  
    
    input$recalc
    N = isolate(input$n_sim)
    Unit = input$unit
    ## Combine data and run the simulation
    if (Unit == '%ID') {
    if (N == 1) {    
      out <- mod%>%param(c(pars,ppars))%>%
        update(atol = 1e-8, rol=1e-4, maxstep = 5000)%>%
        mrgsim_d (data = ex, tgrid = tsamp)
      
      out.sum = cbind.data.frame(Time    = out$time, 
                                 Blood   = ((out$Blood)/Total_Dose)*100,
                                 Lung    = ((out$Lung)/Total_Dose)*100,
                                 GI      = ((out$GI)/Total_Dose)*100, 
                                 Spleen  = ((out$Spleen)/Total_Dose)*100,
                                 Rest    = ((out$Rest)/Total_Dose)*100,
                                 Liver   = ((out$Liver)/Total_Dose)*100,
                                 Kidney  = ((out$Kidney)/Total_Dose)*100,
                                 Urine   = ((out$Urine)/Total_Dose)*100,
                                 Feces   = ((out$Feces)/Total_Dose)*100)
      
      out.sum <- out.sum %>% filter(Time != 0)
      
    } else {
      set.seed (123) 
      out<- mod %>% param(pars) %>%
        data_set(ex)%>%
        idata_set(idata)%>%
        update(atol = 1e-8, rol=1e-4, maxstep = 5000)%>%
        mrgsim(obsonly=TRUE, tgrid=tsamp)%>%filter(time !=0)
      
      
      out.sum = out %>% as_data_frame %>%
        mutate(Time = time) %>%
        dplyr::select(ID, Time, Blood:Feces) %>%
        arrange(Time)  %>%
        group_by(Time) %>%
        dplyr::summarise(Blood.5   = (quantile(Blood, 0.025)/Total_Dose)*100,
                         Blood.25  = (quantile(Blood, 0.25)/Total_Dose)*100,  
                         Blood.50  = (quantile(Blood, 0.50)/Total_Dose)*100,  
                         Blood.75  = (quantile(Blood, 0.75)/Total_Dose)*100,
                         Blood.95  = (quantile(Blood, 0.975)/Total_Dose)*100,
                         Lung.5    = (quantile(Lung, 0.025)/Total_Dose)*100,
                         Lung.25   = (quantile(Lung, 0.25)/Total_Dose)*100,
                         Lung.50   = (quantile(Lung, 0.50)/Total_Dose)*100,    
                         Lung.75   = (quantile(Lung, 0.75)/Total_Dose)*100,
                         Lung.95   = (quantile(Lung, 0.975)/Total_Dose)*100,
                         GI.5      = (quantile(GI, 0.025)/Total_Dose)*100,     
                         GI.25     = (quantile(GI, 0.25)/Total_Dose)*100,     
                         GI.50     = (quantile(GI, 0.50)/Total_Dose)*100,
                         GI.75     = (quantile(GI, 0.75)/Total_Dose)*100,      
                         GI.95     = (quantile(GI, 0.975)/Total_Dose)*100,
                         Spleen.5  = (quantile(Spleen, 0.025)/Total_Dose)*100, 
                         Spleen.25 = (quantile(Spleen, 0.25)/Total_Dose)*100,  
                         Spleen.50 = (quantile(Spleen, 0.5)/Total_Dose)*100,  
                         Spleen.75 = (quantile(Spleen, 0.75)/Total_Dose)*100,
                         Spleen.95 = (quantile(Spleen, 0.975)/Total_Dose)*100,
                         Rest.5    = (quantile(Rest, 0.025)/Total_Dose)*100,
                         Rest.25   = (quantile(Rest, 0.25)/Total_Dose)*100,   
                         Rest.50   = (quantile(Rest, 0.5)/Total_Dose)*100,
                         Rest.75   = (quantile(Rest, 0.75)/Total_Dose)*100,
                         Rest.95   = (quantile(Rest, 0.975)/Total_Dose)*100,
                         Liver.5   = (quantile(Liver, 0.025)/Total_Dose)*100,  
                         Liver.25  = (quantile(Liver, 0.25)/Total_Dose)*100,  
                         Liver.50  = (quantile(Liver, 0.5)/Total_Dose)*100,     
                         Liver.75  = (quantile(Liver, 0.75)/Total_Dose)*100,     
                         Liver.95  = (quantile(Liver, 0.975)/Total_Dose)*100, 
                         Kidney.5  = (quantile(Kidney, 0.025)/Total_Dose)*100, 
                         Kidney.25 = (quantile(Kidney, 0.25)/Total_Dose)*100,    
                         Kidney.50 = (quantile(Kidney, 0.50)/Total_Dose)*100,    
                         Kidney.75 = (quantile(Kidney, 0.75)/Total_Dose)*100,
                         Kidney.95 = (quantile(Kidney, 0.975)/Total_Dose)*100,
                         Urine.5   = (quantile(Urine, 0.025)/Total_Dose)*100,  
                         Urine.25  = (quantile(Urine, 0.25)/Total_Dose)*100,    
                         Urine.50  = (quantile(Urine, 0.5)/Total_Dose)*100,    
                         Urine.75  = (quantile(Urine, 0.75)/Total_Dose)*100,    
                         Urine.95  = (quantile(Urine, 0.975)/Total_Dose)*100,
                         Feces.5   = (quantile(Feces, 0.025)/Total_Dose)*100,  
                         Feces.25  = (quantile(Feces, 0.025)/Total_Dose)*100,  
                         Feces.50  = (quantile(Feces, 0.5)/Total_Dose)*100,     
                         Feces.75  = (quantile(Feces, 0.75)/Total_Dose)*100,
                         Feces.95  = (quantile(Feces, 0.975)/Total_Dose)*100)
      
      
        }
    } else {
      VB = 0.0717*BW
      VLu = 0.0074*BW
      VGI = 0.104*BW
      VS  = 0.0038*BW
      VL  = 0.0346*BW
      VK  = 0.0091*BW
      VRet = BW*(1-0.0717-0.0074-0.104-0.0038-0.0346-0.0091)

      if (N == 1) {    
        out <- mod%>%param(c(pars,ppars))%>%
          update(atol = 1e-8, rol=1e-4, maxstep = 5000)%>%
          mrgsim_d (data = ex, tgrid = tsamp)
        out.sum = cbind.data.frame(Time    = out$time, 
                                   Blood   = (out$Blood)/VB,
                                   Lung    = (out$Lung)/VLu,
                                   GI      = (out$GI)/VGI, 
                                   Spleen  = (out$Spleen)/VS,
                                   Rest    = (out$Rest)/VRet,
                                   Liver   = (out$Liver)/VL,
                                   Kidney  = (out$Kidney)/VK)

        out.sum <- out.sum %>% filter(Time != 0)
        
      } else {
        set.seed (123) 
        out<- mod %>% param(pars) %>%
          data_set(ex)%>%
          idata_set(idata)%>%
          update(atol = 1e-8, rol=1e-4, maxstep = 5000)%>%
          mrgsim(obsonly=TRUE, tgrid=tsamp)%>%filter(time !=0)
        
        
        out.sum = out %>% as_data_frame %>%
          mutate(Time = time) %>%
          dplyr::select(ID, Time, Blood:Feces) %>%
          arrange(Time)  %>%
          group_by(Time) %>%
          dplyr::summarise(Blood.5   = (quantile(Blood, 0.025)/VB),
                           Blood.25  = (quantile(Blood, 0.25)/VB),  
                           Blood.50  = (quantile(Blood, 0.50)/VB),  
                           Blood.75  = (quantile(Blood, 0.75)/VB),
                           Blood.95  = (quantile(Blood, 0.975)/VB),
                           Lung.5    = (quantile(Lung, 0.025)/VLu),
                           Lung.25   = (quantile(Lung, 0.25)/VLu),
                           Lung.50   = (quantile(Lung, 0.50)/VLu),    
                           Lung.75   = (quantile(Lung, 0.75)/VLu),
                           Lung.95   = (quantile(Lung, 0.975)/VLu),
                           GI.5      = (quantile(GI, 0.025)/VGI),     
                           GI.25     = (quantile(GI, 0.25)/VGI),     
                           GI.50     = (quantile(GI, 0.50)/VGI),
                           GI.75     = (quantile(GI, 0.75)/VGI),      
                           GI.95     = (quantile(GI, 0.975)/VGI),
                           Spleen.5  = (quantile(Spleen, 0.025)/VS), 
                           Spleen.25 = (quantile(Spleen, 0.25)/VS),  
                           Spleen.50 = (quantile(Spleen, 0.5)/VS),  
                           Spleen.75 = (quantile(Spleen, 0.75)/VS),
                           Spleen.95 = (quantile(Spleen, 0.975)/VS),
                           Rest.5    = (quantile(Rest, 0.025)/VRet),
                           Rest.25   = (quantile(Rest, 0.25)/VRet),   
                           Rest.50   = (quantile(Rest, 0.5)/VRet),
                           Rest.75   = (quantile(Rest, 0.75)/VRet),
                           Rest.95   = (quantile(Rest, 0.975)/VRet),
                           Live.r5   = (quantile(Liver, 0.025)/VL),  
                           Live.r25  = (quantile(Liver, 0.25)/VL),  
                           Live.r50  = (quantile(Liver, 0.5)/VL),     
                           Live.r75  = (quantile(Liver, 0.75)/VL),     
                           Live.r95  = (quantile(Liver, 0.975)/VL), 
                           Kidn.ey5  = (quantile(Kidney, 0.025)/VK), 
                           Kidney.25 = (quantile(Kidney, 0.25)/VK),    
                           Kidney.50 = (quantile(Kidney, 0.50)/VK),    
                           Kidney.75 = (quantile(Kidney, 0.75)/VK),
                           Kidney.95 = (quantile(Kidney, 0.975)/VK))
                        
        
      }
    }
   
    return(out.sum)
  })
  
  plot <- reactive({
    input$action|input$recalc|input$action1
    target = isolate(input$target)
    N = isolate(input$n_sim)
    
    if(input$unit == '%ID') {y_label = "Percentage of injected dose" 
                             Unit = '(%ID)'
    }else {y_label = "Concentrations of AuNPs "
           Unit = '(mg/L)'}
      
    if (N == 1) {
      a = ggplot(r1(), aes_string("Time", target)) + 
        geom_line(color = 'steelblue4', size = 1, linetype = 'solid', show.legend = F) +
        labs(x = "\nTime (hours)", y = paste(y_label, "in",target, Unit, "\n")) +
        theme_bw() + 
        theme(axis.text=element_text(size=16, face = "bold"),
              axis.title = element_text(size=16, face="bold"),
              legend.title = element_blank()) 
    }else {
      a = ggplot(r1(), aes(Time)) + 
        geom_ribbon(aes_string(ymin = paste0(target,".5"), ymax = paste0(target,".95")), fill = 'steelblue',
                    show.legend = T, size = 0.2, alpha = 0.1) +
        geom_ribbon(aes_string(ymin = paste0(target,".25"), ymax = paste0(target,".75")), fill = 'steelblue4',
                    show.legend = T, size = 0.2, alpha = 0.3) +
        geom_line(aes_string(y = paste0(target,".50")), color="steelblue4", size = 1, show.legend = T) + 
        labs(x = "\nTime (hours)", y = paste(y_label, "in",target, Unit, "\n")) +
        theme_bw() + 
        theme(axis.text=element_text(size=16, face = "bold"),
              axis.title = element_text(size=16, face="bold"),
              legend.title = element_blank()) 
    } 
    
    return(a)
    
  })
  
  
  output$tarplot <- renderPlot({
    input$action
    target = isolate(input$target)
    withProgress(message = 'Creating plot', 
                 detail = 'Please Wait..', value = 0.5, 
                 expr = (print(plot()))) 
    
  })
  
  
  #/////////////////////////////////////////
  # Make a table
  
  r2 <- reactive({
    
    # Refine the dose regiment
    input$action
    tinterval = isolate (input$tinterval)
    end_time  = isolate (input$simu_time)
    Tdoses    = isolate (input$numdose)
    
    dt        = 0.1      
    startime  = 0.1
    stoptime  = round(startime + tinterval*(Tdoses-1) + end_time*24)
    nsim      = (stoptime - startime)/ dt + 1
    Timesim   = seq (startime, stoptime, dt)
    dtout     = 1 ## output dt
    oti       = round(seq (startime, nsim, by = dtout/dt)) # output time interval
    
    data <- r1()
    datatable <- as.data.frame(data)
    return(datatable)
  })
  
  output$table1<- renderTable ({
    N = isolate(input$n_sim)
    if (N == 1) {
    datatable <- r2() 
    } else {datatable <- r2() %>% select(ends_with(c('me',".5",".50",".95")))}
    datatable <-datatable %>% select("Time",order(colnames(datatable[,-1])))
  }, width = '100%',digits = -2)
  
  output$downloadtable <- downloadHandler(
    filename = 'Table.csv',
    content = function(file) {
      write.csv(format(as.data.frame(r1()), digits = 4),file)
    }
  )
  
  #///////////////////////////////////////////
  # make the download link
  output$downb <- downloadHandler(
    filename = paste0("Concentrations in", " ", input$target, ".tiff"),
    content = function(file) {
      ggsave(file, plot(), width = 8, height = 4, dpi = 320)
    },
    contentType = 'image/tiff'
  )
  
  
  
}



