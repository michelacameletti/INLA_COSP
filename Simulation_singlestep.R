Simulation_singlestep = function(x, scenario.settings, 
                                name_file_performance, 
                                sigma2_xi,
                                b0, b1, sigma2_e,
                                sigma2iid_Poisson,
                                gamma0_Poisson){

      
      scenario_n = x #number of the scenario (useful for exporting the output)
      
      #--- Set some parameters according to the scenario
      sample.by = scenario.settings[x,"sample.by"] #how to sample monitoring stations
      prop.points.for.area = scenario.settings[x,"prop.points.for.area"]  #proportion of monitoring stats for each area to be sampled
      
      range.scale = scenario.settings[x,"range.scale"]
      res_pred_grid = scenario.settings[x,"res_pred_grid"]
      
      fixed_E = scenario.settings[x,"E"]
        RR = scenario.settings[x,"RR"]
        gamma1_Poisson = log(RR) 
        gamma0_Poisson = scenario.settings[x,"gamma0"]
      
      
      #--- Simulate the covariate values at the regular grid points
      set.seed(123)
      z_sim_grid = rnorm(highres_simulation_reggrid^2)
      
      
      #values of the fixed parameter for the PRIOR exposure method
      if(scenario_n %% 2 > 0) { #odd scenario with RR=1.05
        b <- seq(-0.2, 0.2, by = 0.01 )  # 0.001  --> 401 values.005
      } else { #even scenario with RR=1.25
        b <- seq(-0.5, 0.5, by = 0.01)   #0.001 --> 401 values with 0.0025 #0.001 with 1001 values
      }
      
      
      # ------------------------------------------------------------
      # ---------- FOR LOOP OVER simulation (1,...,n_simulations)
      # ------------------------------------------------------------
      
      posterior_b0_PI_FF = c()
      posterior_b1_PI_FF = c()
      posterior_sigma2e_PI_FF = c()
      posterior_sigma2xi_PI_FF = c()
      posterior_range_PI_FF =  c()
      cor_pred_obs_area_method1 =  cor_pred_obs_area_method2 = c()
      cor_pred_obs_predgrid_method1 =  cor_pred_obs_predgrid_method2 = c()
      posterior_gamma0_PI_method1 = posterior_gamma0_PI_method2 = c()
      posterior_gamma0_FF_method1 = posterior_gamma0_FF_method2 = c() 
      posterior_gamma0_PE_method1 = posterior_gamma0_PE_method2 = c()
      posterior_gamma1_PI_method1 = posterior_gamma1_PI_method2 = c()
      posterior_gamma1_FF_method1 = posterior_gamma1_FF_method2 = c() 
      posterior_gamma1_PE_method1 = posterior_gamma1_PE_method2 = c()
      posterior_sigma2iid_PI_method1 = posterior_sigma2iid_PI_method2 = c()
      posterior_sigma2iid_FF_method1 = posterior_sigma2iid_FF_method2 = c() 
      posterior_sigma2iid_PE_method1 = posterior_sigma2iid_PE_method2 = c()
      
      n.skip = 0
      which.skip = c()
      which.skip.withnames = c()
      
      
      for(sim_i in 1:n_simulations){
            ############################################################
            # Simulate the data
            ############################################################
            simulated.data = Data_simulation(scenario_n,
                                             sample.by,
                                             prop.points.for.area,
                                             res_pred_grid,
                                             range.scale, sigma2_xi, b0, b1, sigma2_e,
                                             z_sim_grid,
                                             gamma1_Poisson, gamma0_Poisson,
                                             sigma2iid_Poisson, fixed_E)
  
            load(paste0(getwd(), 
                        "/Data","/weights_grid_areas_intersection_res_pred_grid_",res_pred_grid,".Rdata"))
            load(paste0(getwd(), 
                        "/Data","/weights_points_areas_inside_closest_",res_pred_grid,".Rdata")) 
            
            #--- Extract elements from the list
            shapefile_sim_complete = simulated.data$shapefile_sim_complete
            n_areas = nrow(shapefile_sim_complete@data)
            
            cat(paste("@@@@@@@@@@ Scenario n.", scenario_n, "-------------------> Simulation n.", sim_i, "of", n_simulations,"\n"))
            cat(paste("N. monitoring stats = ", simulated.data$n.monitoringstats,"\n"))
            
            
            monitoringstations.data = simulated.data$monitoringstations.data
            monitoringstations.coords = simulated.data$monitoringstations.coords
            #NB: monitoring stats coordinates do not change across iteratios (given a scenario)
            sim_grid_inside_df = simulated.data$sim_grid_inside_df 
            sim_grid_inside_df_pred = simulated.data$sim_grid_inside_df_pred 
            
            shapefile_sim_centroids = gCentroid(shapefile_sim_complete, byid = T) 
            
            ##### Uncomment these lines if you want to change the station locations for each iterations
            #file.remove(paste0(getwd(),
            #                    "/Simulations/Data",
            #                    "/monitoringstationsdata_scenario",scenario_n,".Rdata"))
            
            pred_grid = simulated.data$pred_grid
            z_pred_grid = simulated.data$sim_grid_inside_df_pred$z
            
            max_dist = simulated.data$max_dist
            r = max_dist * range.scale #spatial range
            kappa_xi = sqrt(8) / r
            
            if(sim_i==1){
              posterior_lincombs_method1 = vector("list",length=n_simulations)
              posterior_lincombs_method2 = vector("list",length=n_simulations)
              simulated_data_shapefile = vector("list",length=n_simulations)
            }
            
            ############################################################
            # Inla mesh + spatial process definition
            ############################################################
            
            domain = rbind(c(min(shapefile_sim_centroids$x),min(shapefile_sim_centroids$y)),
                           c(max(shapefile_sim_centroids$x),min(shapefile_sim_centroids$y)),
                           c(min(shapefile_sim_centroids$x),max(shapefile_sim_centroids$y)),
                           c(max(shapefile_sim_centroids$x),max(shapefile_sim_centroids$y)))
            
            mesh_sim = inla.mesh.2d(#loc = pred_grid,
              loc.domain=domain,
              max.edge = c(0.15,0.5),
              cutoff = c(0.05,2), offset = c(0.2,0.3)) 
    
            #   pdf(paste0(getwd(),"/Figures/scenario_n",scenario_n,"/mesh.pdf"))
            #   plot(mesh_sim,las=1)
            #   points(domain,col=2,pch=19,cex=2)
            #   plot(shapefile_sim_complete,axes=T,add=T,border=4,lwd=2)
            #   dev.off()
               
            spde <- inla.spde2.matern(mesh=mesh_sim, alpha=2)
            spde$n.spde
            s_index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)
            
            ##############################################################
            # ------ DEFINE THE inla.stack OBJECTS 
            ##############################################################
            A_est <- inla.spde.make.A(mesh=mesh_sim, 
                                      loc=as.matrix(monitoringstations.coords))
            stack_est <- inla.stack(data=list(x=monitoringstations.data$exposure_with_meas_err),
                                    A=list(A_est,1),
                                    effects=list(c(s_index, list(intercept=1)),
                                                 list(z=monitoringstations.data$z)), 
                                    tag="est") 
            
            A_pred <- inla.spde.make.A(mesh=mesh_sim,
                                      loc=pred_grid@coords)
            stack_pred <- inla.stack(data=list(x=NA),
                                      A=list(A_pred,1),
                                      effects=list(c(s_index,list(intercept=1)),
                                                   list(z=z_pred_grid)),
                                      tag="pred")
        
            stack <- inla.stack(stack_est, stack_pred)
            
            ##############################################################
            # Define the linear combinations in the INLA style for method 1: intersections
            ##############################################################
            #dim_lp <- nrow(inla.stack.A(stack)) + ncol(inla.stack.A(stack)) 
            dim_lp = nrow(monitoringstations.coords) + nrow(pred_grid@coords)
            
            index_pred <- inla.stack.index(stack,"pred")$data #NB! numbers start from nrow(monitoringstations.coords)+1
          
            #each row is an area, each column is a component of the linear combination
            lincombs.weights_method1 = matrix(0, nrow = n_areas, ncol = length(pred_grid))  
            
            lc_all_method1 = c() 
            for(i in seq(1, length(unique(weights_grid_areas_intersection$AREAKEY)))) { #this corresponds to n_area
              
              #cat("------------- Area n.",i,"----", "\n" )
              pixel_area = weights_grid_areas_intersection[which(weights_grid_areas_intersection$AREAKEY==unique(weights_grid_areas_intersection$AREAKEY)[i]),]
              
              #Extract from the SP_ID string only the number related to the pixel (remove g and other numbers at the end)  
              pixel_area$index = as.numeric(substring(gsub( " .*$", "", pixel_area$pixel ),2))
              pixel_area = pixel_area[order(pixel_area$index),] 
              
              #Fill in the matrix with the weights (sum by row = 1)
              lincombs.weights_method1[i,(index_pred[pixel_area$index] - nrow(monitoringstations.coords))]  = pixel_area$weights  
              
              lc_vec <- rep(NA,times=dim_lp)
              lc_vec[index_pred][pixel_area$index] <- pixel_area$weights
              lc_predictor <- inla.make.lincomb(APredictor = lc_vec)
              names(lc_predictor) <- paste0("lc_area_ID=",pixel_area$AREAKEY[1])
              assign(paste0("lc_area_ID=",pixel_area$AREAKEY[1]), lc_predictor)
              
              lc_all_method1 = c(lc_all_method1, get(paste0("lc_area_ID=",pixel_area$AREAKEY[1])))   
            }
            
            #Check that the sum by row is 1
            if( sum(apply(lincombs.weights_method1,1,sum)) != n_areas) print ("!!!!!!!! Errror in creating the lincomb weight matrix!!!!!!!")
            
            #how many lincombs for each area
            n_lc_all_method1 = unlist(lapply(lc_all_method1,
                                             function(x) length(x[[1]]$APredictor$weight)))
            
            ##############################################################
            # Define the linear combinations in the INLA style for method 2: points inside/closest
            ##############################################################
            #each row is an area, each column is a component of the linear combination
            lincombs.weights_method2 = matrix(0, nrow = n_areas, ncol = length(pred_grid))  
            
            
            lc_all_method2 = c() 
            for(i in seq(1, length(unique(weights_points_areas_inside_closest$AREAKEY)))) { #this corresponds to n_area
              
              #cat("------------- Area n.",i,"----", "\n" )
              pixel_area = weights_points_areas_inside_closest[which(weights_points_areas_inside_closest$AREAKEY==unique(weights_points_areas_inside_closest$AREAKEY)[i]),]
              
              #Fill in the matrix with the weights (sum by row = 1)
              pixel_area$index = pixel_area$pixel
              lincombs.weights_method2[i,(index_pred[pixel_area$index] - nrow(monitoringstations.coords))]  = pixel_area$weights  
              
              lc_vec <- rep(NA,times=dim_lp)
              lc_vec[index_pred][pixel_area$index] <- pixel_area$weights
              lc_predictor <- inla.make.lincomb(APredictor = lc_vec)
              names(lc_predictor) <- paste0("lc_area_ID=",pixel_area$AREAKEY[1])
              assign(paste0("lc_area_ID=",pixel_area$AREAKEY[1]), lc_predictor)
              
              lc_all_method2 = c(lc_all_method2, get(paste0("lc_area_ID=",pixel_area$AREAKEY[1])))   
            }
            
            #Check that the sum by row is 1
            if( sum(apply(lincombs.weights_method2,1,sum)) != n_areas) print ("!!!!!!!! Errror in creating the lincomb weight matrix!!!!!!!")
            
            #How many lincombs for each area
            n_lc_all_method2 = unlist(lapply(lc_all_method2,function(x) length(x[[1]]$APredictor$weight)))
            
            ##############################################################
            # FIRST STAGE: run SPDE with linear combinations
            ##############################################################
            formula <- x ~ -1 + intercept + z + f(spatial.field, model=spde)
            
            #The two runs differ in the lincomb specification
            start.exposure = Sys.time()
            output_withlincomb_method1 = inla(formula,
                                      data=inla.stack.data(stack, spde=spde),
                                      family="gaussian",
                                      #control.family = list(hyper=list(prec=list(prior = "logtnormal",
                                      #                                 param = c(sqrt(sigma2_e),100)))), #prior on the sd
                                      control.compute=list(config=TRUE),
                                      control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
                                      lincomb = lc_all_method1, verbose=F)
                                      
            stop.exposure = Sys.time()
            
            output_withlincomb_method2 = inla(formula,
                                              data=inla.stack.data(stack, spde=spde),
                                              family="gaussian",
                                              #control.family = list(hyper=list(prec=list(prior = "logtnormal",
                                              #                                 param = c(sqrt(sigma2_e),100)))), #prior on the sd
                                              control.compute=list(config=TRUE),
                                              control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
                                              lincomb = lc_all_method2, verbose=F)
            
            
          output_withlincomb = output_withlincomb_method1 #parameter posterior distributions do not change with method1/method2
          
          #--- Sample n_random_fromposterior values from the parameter posterior distributions (for the exposure model)
          #--- and save the results
          posterior_b0_PI_FF = c(posterior_b0_PI_FF,
                                 inla.rmarginal(n_random_fromposterior, output_withlincomb$marginals.fixed$intercept))
          posterior_b1_PI_FF = c(posterior_b1_PI_FF,
                                 inla.rmarginal(n_random_fromposterior, output_withlincomb$marginals.fixed$z))
          
          meas_err_var = try(inla.tmarginal(function(x) 1/x,
                                        output_withlincomb$marginals.hyper[[1]]))
          if (class(meas_err_var) == "try-error") next
          
          posterior_sigma2e_PI_FF = c(posterior_sigma2e_PI_FF,
                                      inla.rmarginal(n_random_fromposterior, meas_err_var))
          
          output.field <- inla.spde2.result(inla=output_withlincomb, name="spatial.field",
                                            spde=spde, do.transf=TRUE)
          posterior_range_PI_FF = c(posterior_range_PI_FF,
                                    inla.rmarginal(n_random_fromposterior, output.field$marginals.range.nominal[[1]]))
          posterior_sigma2xi_PI_FF = c(posterior_sigma2xi_PI_FF,
                                       inla.rmarginal(n_random_fromposterior, output.field$marginals.variance.nominal[[1]]))
          
          #--- Summaries from the posterior distributions of the lincombs:
          lincomb_marginals_method1 = output_withlincomb_method1$marginals.lincomb.derived
          lincomb_marginals_method2 = output_withlincomb_method2$marginals.lincomb.derived 
          
          
          # # check the prediction at the grid level
          cor_pred_obs_predgrid_method1[sim_i] = cor(sim_grid_inside_df_pred$exposure,
                                                     output_withlincomb_method1$summary.linear.predictor$`0.5quant`[index_pred])
          cor_pred_obs_predgrid_method2[sim_i] = cor(sim_grid_inside_df_pred$exposure,
                                                     output_withlincomb_method2$summary.linear.predictor$`0.5quant`[index_pred])


          # # The linear predictor is the same with method1/method2
          # plot(output_withlincomb_method1$summary.linear.predictor$`0.5quant`[index_pred],
          # output_withlincomb_method2$summary.linear.predictor$`0.5quant`[index_pred])
          
          # This is just a check for the prediction (done only for the last simulation)
          if(sim_i==n_simulations){ 
            
            check_df = simulated.data$sim_grid_inside_df_pred
            check_df$exposure_est = output_withlincomb_method2$summary.linear.predictor[index_pred, "mean"]
            true_area_exposure = aggregate(exposure ~ AREAKEY,
                                           data=check_df,
                                           FUN=mean)
            est_area_exposure = aggregate(exposure_est ~ AREAKEY,
                                           data=check_df,
                                           FUN=mean)
            shapefile_sim_check = merge(shapefile_sim_complete,
                                  true_area_exposure,
                                  by = "AREAKEY",
                                  all.x = T)
            shapefile_sim_check = merge(shapefile_sim_check,
                                        est_area_exposure,
                                        by = "AREAKEY",
                                        all.x = T)
            pdf(paste0(getwd(), "/Figures/Exposurebias_area_scenario",scenario_n,".pdf"),width=5, height=6)  
            p = spplot(shapefile_sim_check,c("true_area_exposure","exposure","exposure_est"))
            print(p)
            #true_area_exposure is different from exposure because the first is computed
            #on the big number of points used for exposure simulation while the latter
            #is compute using the prediciton grid proints
            dev.off()
            
            
            pdf(paste0(getwd(), "/Figures/Exposurebias_predgrid_scenario",scenario_n,".pdf"),width=5, height=6)  
            mean_pred_true = matrix(simulated.data$sim_grid_inside_df_pred$exposure, res_pred_grid, res_pred_grid, byrow=T)
            mean_pred_grid = matrix(output_withlincomb_method2$summary.linear.predictor[index_pred, "mean"], res_pred_grid, res_pred_grid, byrow=T)
            #image.plot(mean_pred_true)
            #image.plot(mean_pred_grid)
          
            mypalette = colorRampPalette(brewer.pal(11,"RdBu"))(100)
            diff_predgrid = mean_pred_grid - mean_pred_true
            image.plot(simulated.data$pred_grid_x1,simulated.data$pred_grid_x2,diff_predgrid,
                     col= mypalette, main=paste ("Scenario=",scenario_n,"Iteration=",sim_i))
            plot(shapefile_sim_complete,add=T,lwd=2)  
            dev.off()
          }
          
          ##############################################################
          # AREA EXPOSURE ESTIMATION METHOD 1: LINCOMBS with intersections
          # Extract EXPOSURE relevant results for PLUG-IN PI (posterior median)
          ##############################################################
          #--- Median + Create data frame with the posterior median + Add to the shapefile
          posterior_median_method1 = unlist(lapply(lincomb_marginals_method1, function(x) inla.qmarginal(0.5, x)))
          
          posterior_median_method1 = data.frame("AREAKEY"=substring(names(posterior_median_method1),
                                                               nchar(names(posterior_median_method1[1]))-5,max(nchar(names(posterior_median_method1)))),
                                           posterior_median_method1)
  
          shapefile_sim_complete = merge(shapefile_sim_complete,
                                              posterior_median_method1,by="AREAKEY",all.x=T)
  
          #spplot(shapefile_sim_complete,c("true_area_exposure","posterior_median_method1"))
          #plot(shapefile_sim_complete@data$true_area_exposure,shapefile_sim_complete@data$posterior_median_method1)
          cor_pred_obs_area_method1[sim_i] = cor(shapefile_sim_complete@data$true_area_exposure,
                                            shapefile_sim_complete@data$posterior_median_method1)
          
          ##############################################################
          # AREA EXPOSURE ESTIMATION METHOD 2: AREA AVERAGE (with points inside/nearest point)
          ##############################################################
          #--- Median + Create data frame with the posterior median + Add to the shapefile
          posterior_median_method2 = unlist(lapply(lincomb_marginals_method2, function(x) inla.qmarginal(0.5, x)))
          
          posterior_median_method2 = data.frame("AREAKEY"=substring(names(posterior_median_method2),
                                                                    nchar(names(posterior_median_method2[1]))-5,max(nchar(names(posterior_median_method2)))),
                                                posterior_median_method2)
          
          shapefile_sim_complete = merge(shapefile_sim_complete,
                                              posterior_median_method2,by="AREAKEY",all.x=T)
          
          
          #spplot(shapefile_sim_complete,c("true_area_exposure","posterior_median_method1", "posterior_median_method2"))
          #plot(shapefile_sim_complete@data$true_area_exposure,shapefile_sim_complete@data$posterior_median_method1)
          #points(shapefile_sim_complete@data$true_area_exposure,shapefile_sim_complete@data$posterior_median_method2,col=2)
          cor_pred_obs_area_method2[sim_i] = cor(shapefile_sim_complete@data$true_area_exposure,shapefile_sim_complete@data$posterior_median_method2)
          
          
          ###############################################################
          ## Simulate from EXPOSURE posterior distribution at the area level: method1
          ###############################################################
          # Sample n_random_fromPPD samples from the joint posterior and compute manually the lincombs
          # (These samples will be used for computing performance indexes)
          
          list_samples_joint_posterior_method1 = list()
          list_samples_joint_posterior_method2 = list()
          
          for(jj in 1:n_random_fromPPD){
            sample_joint_posterior_method1 = inla.posterior.sample(1,output_withlincomb_method1)
            sample_joint_posterior_method2 = inla.posterior.sample(1,output_withlincomb_method2)
            
            # Extract the part about the predictions (at the regular grid resolution)
            sample_joint_posterior_grid_method1 = sample_joint_posterior_method1[[1]]$latent[1:dim_lp][index_pred]
            sample_joint_posterior_grid_method2 = sample_joint_posterior_method2[[1]]$latent[1:dim_lp][index_pred]
            #length(sample_joint_posterior_grid_method1) #this has the same resolution of the regular grid!
            
            # Compute the area prediction for method1 
            sample_joint_posterior_area_method1 = c()
            for(j in 1:n_areas){
              index_area = lc_all_method1[[j]][[1]]$APredictor$idx - nrow(monitoringstations.coords)
              #NB: indexes start from  nrow(monitoringstations.coords)+1
              #that's the reason why we have - nrow(monitoringstations.coords)
              weight_area = lc_all_method1[[j]][[1]]$APredictor$weight
              sample_joint_posterior_area_method1[j] = sum(sample_joint_posterior_grid_method1[index_area] * weight_area) #posterior area exposure mean
            }
            
            
            list_samples_joint_posterior_method1[[jj]] = data.frame(sample_joint_post_exposure_method1 = sample_joint_posterior_area_method1,
                                                             "AREAKEY"=unique(weights_grid_areas_intersection$AREAKEY))
            
            sample_joint_posterior_area_method2 = c()
            for(j in 1:n_areas){
              index_area = lc_all_method2[[j]][[1]]$APredictor$idx - nrow(monitoringstations.coords)
              #NB: indexes start from  nrow(monitoringstations.coords)+1
              #that's the reason why we have - nrow(monitoringstations.coords)
              weight_area = lc_all_method2[[j]][[1]]$APredictor$weight
              sample_joint_posterior_area_method2[j] = sum(sample_joint_posterior_grid_method2[index_area] * weight_area)
            }
            
            
            
            list_samples_joint_posterior_method2[[jj]] = data.frame(sample_joint_post_exposure_method2 = sample_joint_posterior_area_method2,
                                                                     "AREAKEY"=unique(weights_points_areas_inside_closest$AREAKEY))
             
          }
        

          #Save the same simulated values in another format for computing at the end the performance indexes
          temp_method1 = as.data.frame(lapply(list_samples_joint_posterior_method1,
                                               function(x) x$sample_joint_post_exposure))
          colnames(temp_method1) = paste("Sim",seq(1,n_random_fromPPD))
          temp_method1 = as.matrix(temp_method1)
          posterior_lincombs_method1[[sim_i]] = temp_method1 #matrix n_areas * n_random_fromPPD
          
          #Save the same simulated values in another format for computing at the end the performance indexes
          temp_method2 = as.data.frame(lapply(list_samples_joint_posterior_method2, function(x) x$sample_joint_post_exposure))
          colnames(temp_method2) = paste("Sim",seq(1,n_random_fromPPD))
          temp_method2 = as.matrix(temp_method2)
          posterior_lincombs_method2[[sim_i]] = temp_method2 #matrix n_areas * n_random_fromPPD
          
          
          
          ##############################################################################################################################
          ## SECOND STAGE POISSON MODEL with PI approach (method1 and method2)
          ##############################################################################################################################
          #--- Prepare the data
          data_Poisson = shapefile_sim_complete@data
          data_Poisson$id = seq(1:n_areas)
          
          #--- Method1
          formula.Poisson.PI.method1 = y_Poisson ~ 1 + posterior_median_method1 + f(id,model="iid") 
          
          start.Poisson.PI.method1 = Sys.time()
          output.Poisson.PI.method1 = inla(formula.Poisson.PI.method1,family="poisson",
                                    data=data_Poisson,
                                    offset=log(expected))
          
          posterior_gamma0_PI_method1 = c(posterior_gamma0_PI_method1,
                                  inla.rmarginal(n_random_fromposterior, output.Poisson.PI.method1$marginals.fixed$`(Intercept)`))
          posterior_gamma1_PI_method1 = c(posterior_gamma1_PI_method1,
                                  inla.rmarginal(n_random_fromposterior, output.Poisson.PI.method1$marginals.fixed$posterior_median_method1))
          
          sigma2iid_marg_post_PI_method1 = inla.tmarginal(function(x) 1/x,
                                                  output.Poisson.PI.method1$marginals.hyperpar$`Precision for id`)
          posterior_sigma2iid_PI_method1 = c(posterior_sigma2iid_PI_method1,
                                     inla.rmarginal(n_random_fromposterior, sigma2iid_marg_post_PI_method1))
          
          stop.Poisson.PI.method1 = Sys.time()
          
          #--- Method2
          formula.Poisson.PI.method2 = y_Poisson ~ 1 + posterior_median_method2 + f(id,model="iid") 
          
          start.Poisson.PI.method2 = Sys.time()
          output.Poisson.PI.method2 = inla(formula.Poisson.PI.method2,family="poisson",
                                           data=data_Poisson,
                                           offset=log(expected))
          
          posterior_gamma0_PI_method2 = c(posterior_gamma0_PI_method2,
                                          inla.rmarginal(n_random_fromposterior, output.Poisson.PI.method2$marginals.fixed$`(Intercept)`))
          posterior_gamma1_PI_method2 = c(posterior_gamma1_PI_method2,
                                          inla.rmarginal(n_random_fromposterior, output.Poisson.PI.method2$marginals.fixed$posterior_median_method2))
          
          sigma2iid_marg_post_PI_method2 = inla.tmarginal(function(x) 1/x,
                                                          output.Poisson.PI.method2$marginals.hyperpar$`Precision for id`)
          posterior_sigma2iid_PI_method2 = c(posterior_sigma2iid_PI_method1,
                                             inla.rmarginal(n_random_fromposterior, sigma2iid_marg_post_PI_method2))
          
          stop.Poisson.PI.method2 = Sys.time()
          
          
          ##############################################################
          # SECOND STAGE: Poisson model with FF approach (method1 and method2)
          ##############################################################
          #--- Method1
          formula.Poisson.FF.method1 = y_Poisson ~ 1 + sample_joint_post_exposure_method1 + f(id,model="iid") 
          
          start.Poisson.FF.meethod1 = Sys.time()
          for(j_FF in 1:n_random_fromPPD){
            output.Poisson.FF.method1 <- inla(formula.Poisson.FF.method1,family="poisson",
                                      data=data.frame(data_Poisson,
                                                      list_samples_joint_posterior_method1[[j_FF]]),
                                      offset=log(expected))
            
            posterior_gamma0_FF_method1 = c(posterior_gamma0_FF_method1,
                                    inla.rmarginal(n_random_fromposterior, output.Poisson.FF.method1$marginals.fixed$`(Intercept)`))
            posterior_gamma1_FF_method1 = c(posterior_gamma1_FF_method1,
                                    inla.rmarginal(n_random_fromposterior, output.Poisson.FF.method1$marginals.fixed$sample_joint_post_exposure_method1))
            
            sigma2iid_marg_post_FF.method1 = inla.tmarginal(function(x) 1/x,
                                                    output.Poisson.FF.method1$marginals.hyperpar$`Precision for id`)
            posterior_sigma2iid_FF_method1 = c(posterior_sigma2iid_FF_method1,
                                       inla.rmarginal(n_random_fromposterior, sigma2iid_marg_post_FF.method1))
          }
          stop.Poisson.FF.method1 = Sys.time()
          
          #--- Method2
          formula.Poisson.FF.method2 = y_Poisson ~ 1 + sample_joint_post_exposure_method2 + f(id,model="iid") 
          
          start.Poisson.FF.meethod2 = Sys.time()
          for(j_FF in 1:n_random_fromPPD){
            output.Poisson.FF.method2 <- inla(formula.Poisson.FF.method2,family="poisson",
                                              data=data.frame(data_Poisson,
                                                              list_samples_joint_posterior_method2[[j_FF]]),
                                              offset=log(expected))
            
            posterior_gamma0_FF_method2 = c(posterior_gamma0_FF_method2,
                                            inla.rmarginal(n_random_fromposterior, output.Poisson.FF.method2$marginals.fixed$`(Intercept)`))
            posterior_gamma1_FF_method2 = c(posterior_gamma1_FF_method2,
                                            inla.rmarginal(n_random_fromposterior, output.Poisson.FF.method2$marginals.fixed$sample_joint_post_exposure_method2))
            
            sigma2iid_marg_post_FF.method2 = inla.tmarginal(function(x) 1/x,
                                                            output.Poisson.FF.method2$marginals.hyperpar$`Precision for id`)
            posterior_sigma2iid_FF_method2 = c(posterior_sigma2iid_FF_method2,
                                               inla.rmarginal(n_random_fromposterior, sigma2iid_marg_post_FF.method2))
          }
          stop.Poisson.FF.method2 = Sys.time()
          
          ##############################################################
          # SECOND STAGE: Poisson model with PE approach 
          ##############################################################
          
          #--- Define prior mean and precision for exposure for each area (taken from the lincombs marginal posterior distributions)
          all_post_mean_method1 = unlist(lapply(lincomb_marginals_method1, function(x) inla.zmarginal(x,silent=T)$mean))
          all_post_sd_method1 = unlist(lapply(lincomb_marginals_method1, function(x) inla.zmarginal(x,silent=T)$sd))
          all_post_mean_method2 = unlist(lapply(lincomb_marginals_method2, function(x) inla.zmarginal(x,silent=T)$mean))
          all_post_sd_method2 = unlist(lapply(lincomb_marginals_method2, function(x) inla.zmarginal(x,silent=T)$sd))
          
          #data_Poisson$meanexp = apply(cbind(all_post_mean_method1,all_post_mean_method2), 1 ,mean)
          data_Poisson$meanexp_method1 = all_post_mean_method1
          data_Poisson$meanexp_method2 = all_post_mean_method2
          #data_Poisson$precexp = apply(cbind(1/(all_post_sd_method1^2), 1/(all_post_sd_method2^2)), 1 ,mean)
          data_Poisson$precexp_method1 = 1/(all_post_sd_method1^2)
          data_Poisson$precexp_method2 = 1/(all_post_sd_method2^2)
          
          # Another id for purely iid random effect
          data_Poisson$id2 <- data_Poisson$id
          
          # NB: the values of b are set at the beginning of the file (see around line 22)
          
          
          
          
          #--- Method1
          #--- Run INLA in parallel (the number of models is given by length(b))
          start.Poisson.PE.method1 = Sys.time()
          
          PE_results_method1 = mclapply(b, function(X){fit.inla.method1(d=data_Poisson, X)})
          mliks_method1 <- unlist(lapply(PE_results_method1, function(X){X$mlik[1,1]}))
          bma.w.method1 <- exp(mliks_method1 - max(mliks_method1) + dnorm(b, mean = 0, sd = sqrt(1/.001), log = TRUE))
          #Index to remove models with a very low weight
          bma.idx.method1 <- which(bma.w.method1 > 1e-10)
          
          #--- Extract results for gamma1
          gamma1_PE_method1 <- inla.smarginal(cbind(b, bma.w.method1))
          
          #--- Re-scale BMA weights
          bma.w.method1 <- bma.w.method1 /sum(bma.w.method1)
          #--- Re-scale marginal to REALLY integrate 1
          h.method1 <- gamma1_PE_method1$x[2] - gamma1_PE_method1$x[1]
          w.method1 <- sum(h.method1 * gamma1_PE_method1$y)
          gamma1_PE_method1$y <-  gamma1_PE_method1$y / w.method1
          
          #Plot
          #plot(gamma1_PE_method1, type = "l")
          #abline(v=log(RR))
          
          #--- Simulate from the posterior of gamma1
          posterior_gamma1_PE_method1 = c(posterior_gamma1_PE_method1,
                                  inla.rmarginal(n_random_fromposterior, gamma1_PE_method1))
          
          #--- Extract results for intercept
          gamma0_marg_post_PE_method1_list = lapply(PE_results_method1,
                                            function(Y){inla.tmarginal(identity, Y$marginals.fixed$`(Intercept)`)})
          gamma0_PE_method1 <- try(INLABMA:::fitmargBMA(gamma0_marg_post_PE_method1_list[bma.idx.method1], bma.w.method1[bma.idx.method1]))
          if(class(gamma0_PE_method1)=="data.frame"){
            posterior_gamma0_PE_method1 = c(posterior_gamma0_PE_method1,
                                            inla.rmarginal(n_random_fromposterior, gamma0_PE_method1))
          } else { #error in integrate
            n.skip = n.skip + 1
            which.skip.withnames = c(which.skip.withnames,paste0("gamma0_",sim_i))
            which.skip.withnames = c(which.skip.withnames,sim_i)
            
            #Remove for FF, PI, PE and the other parameters the simulated values of this iteration that we are skipping
            posterior_gamma0_FF_method1 = posterior_gamma0_FF_method1[- ((length(posterior_gamma0_FF_method1)-n_random_fromposterior+1):length(posterior_gamma0_FF_method1))]
            posterior_gamma1_FF_method1 = posterior_gamma1_FF_method1[- ((length(posterior_gamma1_FF_method1)-n_random_fromposterior+1):length(posterior_gamma1_FF_method1))]
            posterior_sigma2iid_FF_method1 = posterior_sigma2iid_FF_method1[- ((length(posterior_sigma2iid_FF_method1)-n_random_fromposterior+1):length(posterior_sigma2iid_FF_method1))]
            
            posterior_gamma0_PI_method1 = posterior_gamma0_PI_method1[- ((length(posterior_gamma0_PI_method1)-n_random_fromposterior+1):length(posterior_gamma0_PI_method1))]
            posterior_gamma1_PI_method1 = posterior_gamma1_PI_method1[- ((length(posterior_gamma1_PI_method1)-n_random_fromposterior+1):length(posterior_gamma1_PI_method1))]
            posterior_sigma2iid_PI_method1 = posterior_sigma2iid_PI_method1[- ((length(posterior_sigma2iid_PI_method1)-n_random_fromposterior+1):length(posterior_sigma2iid_PI_method1))]
            
            posterior_gamma1_PE_method1 = posterior_gamma1_PE_method1[- ((length(posterior_gamma1_PE_method1)-n_random_fromposterior+1):length(posterior_gamma1_PE_method1))]
            
            next
          }  
          
          #--- Extract results for sigma2iid
          sigma2iid_marg_post_PE_method1_list = lapply(PE_results_method1,
                                                      function(Y){inla.tmarginal(function(x) 1/x, Y$marginals.hyperpar$`Precision for id2`)})
          
          sigma2iid_PE_method1 = try(INLABMA:::fitmargBMA(sigma2iid_marg_post_PE_method1_list[bma.idx.method1], bma.w.method1[bma.idx.method1]))
          if(class(sigma2iid_PE_method1)=="data.frame"){
            posterior_sigma2iid_PE_method1 = c(posterior_sigma2iid_PE_method1,
                                              inla.rmarginal(n_random_fromposterior, sigma2iid_PE_method1))
          } else { #error in integrate
            n.skip = n.skip + 1
            which.skip.withnames = c(which.skip.withnames,paste0("sigma2iid_",sim_i))
            which.skip = c(which.skip,sim_i)
            
            #Remove for FF, PI, PE and the other parameters the simulated values of this iteration that we are skipping
            posterior_gamma0_FF_method1 = posterior_gamma0_FF_method1[- ((length(posterior_gamma0_FF_method1)-n_random_fromposterior+1):length(posterior_gamma0_FF_method1))]
            posterior_gamma1_FF_method1 = posterior_gamma1_FF_method1[- ((length(posterior_gamma1_FF_method1)-n_random_fromposterior+1):length(posterior_gamma1_FF_method1))]
            posterior_sigma2iid_FF_method1 = posterior_sigma2iid_FF_method1[- ((length(posterior_sigma2iid_FF_method1)-n_random_fromposterior+1):length(posterior_sigma2iid_FF_method1))]
            
            posterior_gamma0_PI_method1 = posterior_gamma0_PI_method1[- ((length(posterior_gamma0_PI_method1)-n_random_fromposterior+1):length(posterior_gamma0_PI_method1))]
            posterior_gamma1_PI_method1 = posterior_gamma1_PI_method1[- ((length(posterior_gamma1_PI_method1)-n_random_fromposterior+1):length(posterior_gamma1_PI_method1))]
            posterior_sigma2iid_PI_method1 = posterior_sigma2iid_PI_method1[- ((length(posterior_sigma2iid_PI_method1)-n_random_fromposterior+1):length(posterior_sigma2iid_PI_method1))]
            
            posterior_gamma0_PE_method1 = posterior_gamma0_PE_method1[- ((length(posterior_gamma0_PE_method1)-n_random_fromposterior+1):length(posterior_gamma0_PE_method1))]
            posterior_gamma1_PE_method1 = posterior_gamma1_PE_method1[- ((length(posterior_gamma1_PE_method1)-n_random_fromposterior+1):length(posterior_gamma1_PE_method1))]

            next
          }
          stop.Poisson.PE.method1 = Sys.time()
          
          
          
          #--- Method2
          #--- Run INLA in parallel (the number of models is given by length(b))
          start.Poisson.PE.method2 = Sys.time()
          PE_results_method2 = mclapply(b, function(X){fit.inla.method2(d=data_Poisson, X)})
          mliks_method2 <- unlist(lapply(PE_results_method2, function(X){X$mlik[1,1]}))
          bma.w.method2 <- exp(mliks_method2 - max(mliks_method2) + dnorm(b, mean = 0, sd = sqrt(1/.001), log = TRUE))
          #Index to remove models with a very low weight
          bma.idx.method2 <- which(bma.w.method2 > 1e-10)
          
          #--- Extract results for gamma1
          gamma1_PE_method2 <- inla.smarginal(cbind(b, bma.w.method2))
          
          #--- Re-scale BMA weights
          bma.w.method2 <- bma.w.method2 /sum(bma.w.method2)
          #--- Re-scale marginal to REALLY integrate 1
          h.method2 <- gamma1_PE_method2$x[2] - gamma1_PE_method2$x[1]
          w.method2 <- sum(h.method2 * gamma1_PE_method2$y)
          gamma1_PE_method2$y <-  gamma1_PE_method2$y / w.method2
          
          #Plot
          #plot(gamma1_PE_method2, type = "l")
          #abline(v=log(RR))
          
          #--- Simulate from the posterior of gamma1
          posterior_gamma1_PE_method2 = c(posterior_gamma1_PE_method2,
                                          inla.rmarginal(n_random_fromposterior, gamma1_PE_method2))
          
          #--- Extract results for intercept
          gamma0_marg_post_PE_method2_list = lapply(PE_results_method2,
                                                    function(Y){inla.tmarginal(identity, Y$marginals.fixed$`(Intercept)`)})
          gamma0_PE_method2 <- try(INLABMA:::fitmargBMA(gamma0_marg_post_PE_method2_list[bma.idx.method2], bma.w.method2[bma.idx.method2]))
          if(class(gamma0_PE_method2)=="data.frame"){
            posterior_gamma0_PE_method2 = c(posterior_gamma0_PE_method2,
                                            inla.rmarginal(n_random_fromposterior, gamma0_PE_method2))
          } else { #error in integrate
            n.skip = n.skip + 1
            which.skip.withnames = c(which.skip.withnames,paste0("gamma0_",sim_i))
            which.skip.withnames = c(which.skip.withnames,sim_i)
            
            #Remove for FF, PI, PE and the other parameters the simulated values of this iteration that we are skipping
            posterior_gamma0_FF_method2 = posterior_gamma0_FF_method2[- ((length(posterior_gamma0_FF_method2)-n_random_fromposterior+1):length(posterior_gamma0_FF_method2))]
            posterior_gamma1_FF_method2 = posterior_gamma1_FF_method2[- ((length(posterior_gamma1_FF_method2)-n_random_fromposterior+1):length(posterior_gamma1_FF_method2))]
            posterior_sigma2iid_FF_method2 = posterior_sigma2iid_FF_method2[- ((length(posterior_sigma2iid_FF_method2)-n_random_fromposterior+1):length(posterior_sigma2iid_FF_method2))]
            
            posterior_gamma0_PI_method2 = posterior_gamma0_PI_method2[- ((length(posterior_gamma0_PI_method2)-n_random_fromposterior+1):length(posterior_gamma0_PI_method2))]
            posterior_gamma1_PI_method2 = posterior_gamma1_PI_method2[- ((length(posterior_gamma1_PI_method2)-n_random_fromposterior+1):length(posterior_gamma1_PI_method2))]
            posterior_sigma2iid_PI_method2 = posterior_sigma2iid_PI_method2[- ((length(posterior_sigma2iid_PI_method2)-n_random_fromposterior+1):length(posterior_sigma2iid_PI_method2))]
            
            posterior_gamma1_PE_method2 = posterior_gamma1_PE_method2[- ((length(posterior_gamma1_PE_method2)-n_random_fromposterior+1):length(posterior_gamma1_PE_method2))]
            next
          }  
          
          #--- Extract results for sigma2iid
          sigma2iid_marg_post_PE_method2_list = lapply(PE_results_method2,
                                                       function(Y){inla.tmarginal(function(x) 1/x, Y$marginals.hyperpar$`Precision for id2`)})
          
          sigma2iid_PE_method2 = try(INLABMA:::fitmargBMA(sigma2iid_marg_post_PE_method2_list[bma.idx.method2], bma.w.method2[bma.idx.method2]))
          if(class(sigma2iid_PE_method2)=="data.frame"){
            posterior_sigma2iid_PE_method2 = c(posterior_sigma2iid_PE_method2,
                                               inla.rmarginal(n_random_fromposterior, sigma2iid_PE_method2))
          } else { #error in integrate
            n.skip = n.skip + 1
            which.skip.withnames = c(which.skip.withnames,paste0("sigma2iid_",sim_i))
            which.skip = c(which.skip,sim_i)
            
            #Remove for FF, PI, PE and the other parameters the simulated values of this iteration that we are skipping
            posterior_gamma0_FF_method2 = posterior_gamma0_FF_method2[- ((length(posterior_gamma0_FF_method2)-n_random_fromposterior+1):length(posterior_gamma0_FF_method2))]
            posterior_gamma1_FF_method2 = posterior_gamma1_FF_method2[- ((length(posterior_gamma1_FF_method2)-n_random_fromposterior+1):length(posterior_gamma1_FF_method2))]
            posterior_sigma2iid_FF_method2 = posterior_sigma2iid_FF_method2[- ((length(posterior_sigma2iid_FF_method2)-n_random_fromposterior+1):length(posterior_sigma2iid_FF_method2))]
            
            posterior_gamma0_PI_method2 = posterior_gamma0_PI_method2[- ((length(posterior_gamma0_PI_method2)-n_random_fromposterior+1):length(posterior_gamma0_PI_method2))]
            posterior_gamma1_PI_method2 = posterior_gamma1_PI_method2[- ((length(posterior_gamma1_PI_method2)-n_random_fromposterior+1):length(posterior_gamma1_PI_method2))]
            posterior_sigma2iid_PI_method2 = posterior_sigma2iid_PI_method2[- ((length(posterior_sigma2iid_PI_method2)-n_random_fromposterior+1):length(posterior_sigma2iid_PI_method2))]
            
            posterior_gamma0_PE_method2 = posterior_gamma0_PE_method2[- ((length(posterior_gamma0_PE_method2)-n_random_fromposterior+1):length(posterior_gamma0_PE_method2))]
            posterior_gamma1_PE_method2 = posterior_gamma1_PE_method2[- ((length(posterior_gamma1_PE_method2)-n_random_fromposterior+1):length(posterior_gamma1_PE_method2))]
            
            next
          }
          stop.Poisson.PE.method2 = Sys.time()
          
          
          ###############################################################
          ## Some final savings
          ###############################################################
          
          
          #--- Save the shapefile data
          simulated_data_shapefile[[sim_i]] = shapefile_sim_complete@data
          
          
          #--- Save a file containing information about each area (how many lincombs and monitoring stats)
          
          n_lc_all = data.frame(n_lc_all_method1,n_lc_all_method2)
          #n_lc_all$AREAKEY = substring(rownames(n_lc_all),
          #              nchar(rownames(n_lc_all[1]))-5,max(nchar(rownames(n_lc_all))))
          n_lc_all$AREAKEY = sub('.*=', '', rownames(n_lc_all)) #extract everything after the =
          
          #--- Add the number of monitoring stations
          nstats = data.frame(table(monitoringstations.data$AREAKEY),
                     AREAKEY=names(table(monitoringstations.data$AREAKEY)))
          nstats = nstats[,-1]
          colnames(nstats)= c("n_monit_stats","AREAKEY")
          info_area = merge(n_lc_all,nstats, all.x=T)
          
          #--- Add the number of points used for simulating the area true level
          nsimpoints = data.frame(table(sim_grid_inside_df$AREAKEY),
                              AREAKEY=names(table(sim_grid_inside_df$AREAKEY)))
          nsimpoints = nsimpoints[,-1]
          colnames(nsimpoints)= c("n_sim_points","AREAKEY")
          info_area2 = merge(info_area,nsimpoints)
          
          
          write.table(info_area2,
                    paste0(getwd(),"/Output/info_area_scenario_",scenario_n,".csv"),
                    sep=";",append=F,row.names = F)
        
      } #close for(sim_i in 1:n_simulations)

    
      ##############################################################
      # Save all the posteriors
      ##############################################################
      correlations = data.frame(rbind(cor_pred_obs_predgrid_method1,
                           cor_pred_obs_predgrid_method2,
                           cor_pred_obs_area_method1,
                           cor_pred_obs_area_method2))
      colnames(correlations) = paste0("sim",seq(1,n_simulations))
      correlations$scenario = rep(scenario_n,4)
      
      
      if(!file.exists(paste0(getwd(),"/Output/correlations_",scenario_n,".csv"))) {
        write.table(correlations,
                  paste0(getwd(),"/Output/correlations_",scenario_n,".csv"),
                  sep=";",append=F,row.names = T)
      } else {
        write.table(correlations,
                    paste0(getwd(),"/Output/correlations_",scenario_n,".csv"),
                    sep=";",append=T,row.names = T,col.names = F)
      }
      
      
      time.exposure = stop.exposure - start.exposure
      
      time.Poisson.PI.method1 = stop.Poisson.PI.method1 - start.Poisson.PI.method1
      time.Poisson.PI.method2 = stop.Poisson.PI.method2 - start.Poisson.PI.method2
      
      time.Poisson.FF.method1 = stop.Poisson.FF.method1 - start.Poisson.FF.meethod1
      time.Poisson.FF.method2 = stop.Poisson.FF.method2 - start.Poisson.FF.meethod2
      
      time.Poisson.PE.method1 = stop.Poisson.PE.method1 - start.Poisson.PE.method1
      time.Poisson.PE.method2 = stop.Poisson.PE.method2 - start.Poisson.PE.method2
      
      save(posterior_gamma0_PI_method1, posterior_gamma0_PI_method2,
           posterior_gamma0_PE_method1, posterior_gamma0_PE_method2,
           posterior_gamma0_FF_method1, posterior_gamma0_FF_method2,
           posterior_gamma1_PI_method1, posterior_gamma1_PI_method2,
           posterior_gamma1_PE_method1, posterior_gamma1_PE_method2,
           posterior_gamma1_FF_method1, posterior_gamma1_FF_method2,
           posterior_sigma2iid_PI_method1,posterior_sigma2iid_PI_method2,
           posterior_sigma2iid_PE_method1, posterior_sigma2iid_PE_method2,
           posterior_sigma2iid_FF_method1,posterior_sigma2iid_FF_method2,
           simulated_data_shapefile,
           #cor_pred_obs_predgrid_method1,cor_pred_obs_predgrid_method2,
           #cor_pred_obs_area_method1,cor_pred_obs_area_method2,
           posterior_lincombs_method1,
           posterior_lincombs_method2,
           sim_grid_inside_df,
           sim_grid_inside_df_pred,
           time.exposure, 
           time.Poisson.PI.method2, time.Poisson.PI.method1,
           time.Poisson.FF.method1, time.Poisson.FF.method2,
           time.Poisson.PE.method1, time.Poisson.PE.method2,
           file=paste0(getwd(),"/Output/","AllposteriorsANDsimulatedData_scenarioN",scenario_n,".RData"))
      
  
  
  ##############################################################
  # Compute bias, RMSE and coverage for the exposure model parameters
  ##############################################################

  #--- Parameter bias + rmse + coverage 
  b0_perf_indexes_PI_FF = unlist(my_bias_RMSE_coverage(posterior_b0_PI_FF, n_simulations, rep(b0,n_simulations)))
  b1_perf_indexes_PI_FF = unlist(my_bias_RMSE_coverage(posterior_b1_PI_FF, n_simulations, rep(b1,n_simulations)))
  sigma2e_perf_indexes_PI_FF = unlist(my_bias_RMSE_coverage(posterior_sigma2e_PI_FF, n_simulations, rep(sigma2_e,n_simulations)))
  sigma2xi_perf_indexes_PI_FF = unlist(my_bias_RMSE_coverage(posterior_sigma2xi_PI_FF, n_simulations, rep(sigma2_xi,n_simulations)))
  range_perf_indexes_PI_FF = unlist(my_bias_RMSE_coverage(posterior_range_PI_FF, n_simulations, rep(r,n_simulations)))
  
  #--- Export exposure model parameter performance
  output_exposurepar_performance = rbind(
    c(b0_perf_indexes_PI_FF,
      b1_perf_indexes_PI_FF,
      sigma2e_perf_indexes_PI_FF,
      sigma2xi_perf_indexes_PI_FF,
      range_perf_indexes_PI_FF)
    )

  output_exposurepar_performance = data.frame(method=c(paste0("PI-FF",scenario_n)),output_exposurepar_performance)
  colnames(output_exposurepar_performance)[-1] = do.call(paste0,
                                                    expand.grid(c("bias","rmse","coverage%"),c("b0","b1","sigma2e","sigma2xi","range")))
  
  if(!file.exists(paste0(getwd(),"/Output/Exposure_",name_file_performance))) {
    write.table(output_exposurepar_performance,
                paste0(getwd(),"/Output/Exposure_",name_file_performance),
                sep=";",row.names = F)
  } else {
    write.table(output_exposurepar_performance,
                paste0(getwd(),"/Output/Exposure_",name_file_performance),
                sep=";",row.names = F,append=T,col.names = F)
  }
  
  ##############################################################
  # Plot posterior distributions EXPOSURE model parameters
  ##############################################################
  pdf(paste0(getwd(),"/Figures/","Posteriors__exposure_scenarioN",scenario_n,".pdf"))
  par(mfrow=c(3,2))
  density_b0_PI_FF = density(posterior_b0_PI_FF)
  plot(density_b0_PI_FF,main="",ylab = expression(b[0]))
  abline(v=b0,lty=4,col=2,lwd=2)
  abline(v=median(posterior_b0_PI_FF))
  boxplot(matrix(posterior_b0_PI_FF,ncol=n_simulations))
  abline(h=b0,lty=4,col=2,lwd=2)
  
  density_b1_PI_FF = density(posterior_b1_PI_FF)
  plot(density_b1_PI_FF,main="",ylab = expression(b[1]))
  abline(v=b1,lty=4,col=2,lwd=2)
  abline(v=median(posterior_b1_PI_FF))
  boxplot(matrix(posterior_b1_PI_FF,ncol=n_simulations))
  abline(h=b1,lty=4,col=2,lwd=2)
  
  density_sigma2e_PI_FF = density(posterior_sigma2e_PI_FF)
  plot(density_sigma2e_PI_FF,main="",ylab = expression(sigma[e]^2))
  abline(v=sigma2_e,lty=4,col=2,lwd=2)
  abline(v=median(posterior_sigma2e_PI_FF))
  boxplot(matrix(posterior_sigma2e_PI_FF,ncol=n_simulations))
  abline(h=sigma2_e,lty=4,col=2,lwd=2)
  
  density_sigma2xi_PI_FF = density(posterior_sigma2xi_PI_FF)
  plot(density_sigma2xi_PI_FF,main="",ylab = expression(sigma[omega]^2))
  abline(v=sigma2_xi,lty=4,col=2,lwd=2)
  abline(v=median(posterior_sigma2xi_PI_FF))
  boxplot(matrix(posterior_sigma2xi_PI_FF,ncol=n_simulations))
  abline(h=sigma2_xi,lty=4,col=2,lwd=2)
  
  density_range_PI_FF = density(posterior_range_PI_FF)
  plot(density_range_PI_FF,main="",ylab = "range")
  abline(v=r,lty=4,col=2,lwd=2)
  abline(v=median(posterior_range_PI_FF))
  boxplot(matrix(posterior_range_PI_FF,ncol=n_simulations))
  abline(h=r,lty=4,col=2,lwd=2)
  
  dev.off()
  
  ##############################################################
  # Compute bias+mse+coverage for exposure area prediction
  ##############################################################

  #comparing methods using the posterior exposure distribution for each area
  list_simulatedvalues_method1 = posterior_lincombs_method1
  list_simulatedvalues_method2 = posterior_lincombs_method2
  
  bias_method1 = bias_method2 = matrix(NA,nrow=n_simulations,ncol=n_areas)
  rmse_method1 = rmse_method2 = matrix(NA,nrow=n_simulations,ncol=n_areas)
  cov_method1 = cov_method2 = matrix(NA,nrow=n_simulations,ncol=n_areas)
  
  for(i in 1:n_simulations){
     temp_bias_method1 = temp_bias_method2 = c()
     temp_cov_method1 = temp_cov_method2 = c()
     temp_rmse_method1 = temp_rmse_method2 = c()
  
     matrix_simulatedvalues_method1 = list_simulatedvalues_method1[[i]]
     matrix_simulatedvalues_method2 = list_simulatedvalues_method2[[i]]
  
     #For each area compute the relative bias and the coverage (nb: the true values change)
     for(A in 1:n_areas){
       #temp_bias[A] = mean((matrix_simulatedvalues[A,] - true_exposure[A]) / true_exposure[A])
       temp_bias_method1[A] = mean(matrix_simulatedvalues_method1[A,] - simulated_data_shapefile[[i]]$true_area_exposure[A])
       temp_bias_method2[A] = mean(matrix_simulatedvalues_method2[A,] - simulated_data_shapefile[[i]]$true_area_exposure[A])
  
       temp_rmse_method1[A] = sqrt(mean((matrix_simulatedvalues_method1[A,] - simulated_data_shapefile[[i]]$true_area_exposure[A])^2))
       temp_rmse_method2[A] = sqrt(mean((matrix_simulatedvalues_method2[A,] - simulated_data_shapefile[[i]]$true_area_exposure[A])^2))
  
       q0025_method1 = quantile(matrix_simulatedvalues_method1[A,],0.025)
       q0975_method1 = quantile(matrix_simulatedvalues_method1[A,],0.975)
       temp_cov_method1[A] = (q0025_method1 < simulated_data_shapefile[[i]]$true_area_exposure[A])  & (simulated_data_shapefile[[i]]$true_area_exposure[A] < q0975_method1)
  
       q0025_method2 = quantile(matrix_simulatedvalues_method2[A,],0.025)
       q0975_method2 = quantile(matrix_simulatedvalues_method2[A,],0.975)
       temp_cov_method2[A] = (q0025_method2 < simulated_data_shapefile[[i]]$true_area_exposure[A])  & (simulated_data_shapefile[[i]]$true_area_exposure[A] < q0975_method2)
     }
  
     bias_method1[i,] = temp_bias_method1
     bias_method2[i,] = temp_bias_method2
  
     cov_method1[i,] = temp_cov_method1
     cov_method2[i,] = temp_cov_method2
  
     rmse_method1[i,] = temp_rmse_method1
     rmse_method2[i,] = temp_rmse_method2
  }
  
   mean_bias_method1 = apply(bias_method1,2,mean)  #mean over simulations
   mean_bias_method2 = apply(bias_method2,2,mean)  #mean over simulations
  
   mean_coverage_method1 = apply(cov_method1,2,sum) / n_simulations * 100
   mean_coverage_method2 = apply(cov_method2,2,sum) / n_simulations * 100
  
   mean_rmse_method1 = apply(rmse_method1,2,mean)
   mean_rmse_method2 = apply(rmse_method2,2,mean)
  
   mean_overiteration_summarystats = data.frame(rbind(mean_bias_method1,mean_bias_method2,
                                                      mean_rmse_method1,mean_rmse_method2,
                                                      mean_coverage_method1,mean_coverage_method2))
   mean_overiteration_summarystats$scenario = rep(scenario_n,6)
   mean_overiteration_summarystats$method = c("mean_bias_method1","mean_bias_method2",
                                              "mean_rmse_method1","mean_rmse_method2",
                                              "mean_coverage_method1","mean_coverage_method2")
  
   if(!file.exists(paste0(getwd(),"/Output/Comparingexposure_meanoveriterations_",name_file_performance))) {
    write.table(mean_overiteration_summarystats,
                 paste0(getwd(),"/Output/Comparingexposure_meanoveriterations_",name_file_performance),
                 sep=";",row.names = F)
   } else {
     write.table(mean_overiteration_summarystats,
                 paste0(getwd(),"/Output/Comparingexposure_meanoveriterations_",name_file_performance),
                 sep=";",row.names = F,append=T,col.names = F)
   }
     
   
   ##############################################################
   # PLOT bias+mse+coverage for exposure area prediction
   ##############################################################
   
  
  #--- Plot the mean bias
  area_bias = data.frame(bias_method1 = mean_bias_method1,bias_method2 = mean_bias_method2)
  area_bias$AREAKEY = shapefile_sim_complete@data$AREAKEY
  shapefile_sim_complete = merge(shapefile_sim_complete, area_bias, by="AREAKEY")

  pdf(paste0(getwd(),"/Figures/","Area_bias_scenarioN",scenario_n,".pdf"))
  mypalette = colorRampPalette(brewer.pal(11,"RdBu"))(100)
  coordinates(monitoringstations.coords) = c("x1","x2")
  proj4string(monitoringstations.coords) = proj4string(shapefile_sim_complete)

  p=spplot(shapefile_sim_complete, c("bias_method1","bias_method2"),
           #names.attr = c("PI-FF"),
           #sp.layout = c('sp.points', monitoringstations.coords, col='magenta', pch=16),
           col.regions=mypalette, cut=16)
           #sp.layout = c('sp.points', pred_grid, col='magenta', pch=16))

  print(p)
  dev.off()

  #--- Plot the mean rmse
  area_RMSE = data.frame(RMSE_method1 = mean_rmse_method1,RMSE_method2 = mean_rmse_method2)
  area_RMSE$AREAKEY = shapefile_sim_complete@data$AREAKEY
  shapefile_sim_complete = merge(shapefile_sim_complete, area_RMSE, by="AREAKEY")
  
  pdf(paste0(getwd(),"/Figures/","Area_rmse_scenarioN",scenario_n,".pdf"))
  mypalette = rev(heat.colors(n=64))

  #coordinates(monitoringstations.coords) = c("x1","x2")
  #proj4string(monitoringstations.coords) = proj4string(shapefile_sim_complete)

  p=spplot(shapefile_sim_complete, c("RMSE_method1","RMSE_method2"),
           col.regions=mypalette, cut=16
           #sp.layout = c('sp.points', monitoringstations.coords, col='magenta', pch=16)
           )
  print(p)
  dev.off()

  #--- Plot the mean coverage
  area_coverage = data.frame(cov_method1 = mean_coverage_method1,cov_method2 = mean_coverage_method2)
  area_coverage$AREAKEY = shapefile_sim_complete@data$AREAKEY
  shapefile_sim_complete = merge(shapefile_sim_complete, area_coverage, by="AREAKEY")
  

  pdf(paste0(getwd(),"/Figures/","Area_perc_coverage_scenarioN",scenario_n,".pdf"))
  mypalette = colorRampPalette(brewer.pal(9,"Blues"))(100)
  p=spplot(shapefile_sim_complete, c("cov_method1","cov_method2"),
           #names.attr = c("PI-FF"),
           col.regions=mypalette, cut=16
           #sp.layout = c('sp.points', monitoringstations.coords, col='magenta', pch=16)
           )

  print(p)
  dev.off()

  
  
  
  ##############################################################
  # Compute bias, RMSE and coverage for the Poisson parameters
  ##############################################################
  #--- Parameter bias + rmse + coverage for the Poisson parameters
  
  gamma0_perf_indexes_PI_method1 = unlist(my_bias_RMSE_coverage(posterior_gamma0_PI_method1, n_simulations, rep(gamma0_Poisson,n_simulations)))
  gamma1_perf_indexes_PI_method1 = unlist(my_bias_RMSE_coverage(posterior_gamma1_PI_method1, n_simulations, rep(gamma1_Poisson,n_simulations)))
  sigma2iid_perf_indexes_PI_method1 = unlist(my_bias_RMSE_coverage(posterior_sigma2iid_PI_method1, n_simulations, rep(sigma2iid_Poisson,n_simulations)))
  
  gamma0_perf_indexes_PI_method2 = unlist(my_bias_RMSE_coverage(posterior_gamma0_PI_method2, n_simulations, rep(gamma0_Poisson,n_simulations)))
  gamma1_perf_indexes_PI_method2 = unlist(my_bias_RMSE_coverage(posterior_gamma1_PI_method2, n_simulations, rep(gamma1_Poisson,n_simulations)))
  sigma2iid_perf_indexes_PI_method2 = unlist(my_bias_RMSE_coverage(posterior_sigma2iid_PI_method2, n_simulations, rep(sigma2iid_Poisson,n_simulations)))
  
  gamma0_perf_indexes_FF_method1 = unlist(my_bias_RMSE_coverage(posterior_gamma0_FF_method1, n_simulations, rep(gamma0_Poisson,n_simulations)))
  gamma1_perf_indexes_FF_method1 = unlist(my_bias_RMSE_coverage(posterior_gamma1_FF_method1, n_simulations, rep(gamma1_Poisson,n_simulations)))
  sigma2iid_perf_indexes_FF_method1 = unlist(my_bias_RMSE_coverage(posterior_sigma2iid_FF_method1, n_simulations, rep(sigma2iid_Poisson,n_simulations)))
  
  gamma0_perf_indexes_FF_method2 = unlist(my_bias_RMSE_coverage(posterior_gamma0_FF_method2, n_simulations, rep(gamma0_Poisson,n_simulations)))
  gamma1_perf_indexes_FF_method2 = unlist(my_bias_RMSE_coverage(posterior_gamma1_FF_method2, n_simulations, rep(gamma1_Poisson,n_simulations)))
  sigma2iid_perf_indexes_FF_method2 = unlist(my_bias_RMSE_coverage(posterior_sigma2iid_FF_method2, n_simulations, rep(sigma2iid_Poisson,n_simulations)))
  
  gamma0_perf_indexes_PE_method1 = unlist(my_bias_RMSE_coverage(posterior_gamma0_PE_method1, n_simulations, rep(gamma0_Poisson,n_simulations)))
  gamma1_perf_indexes_PE_method1 = unlist(my_bias_RMSE_coverage(posterior_gamma1_PE_method1, n_simulations, rep(gamma1_Poisson,n_simulations)))
  sigma2iid_perf_indexes_PE_method1 = unlist(my_bias_RMSE_coverage(posterior_sigma2iid_PE_method1, n_simulations, rep(sigma2iid_Poisson,n_simulations)))
  
  gamma0_perf_indexes_PE_method2 = unlist(my_bias_RMSE_coverage(posterior_gamma0_PE_method2, n_simulations, rep(gamma0_Poisson,n_simulations)))
  gamma1_perf_indexes_PE_method2 = unlist(my_bias_RMSE_coverage(posterior_gamma1_PE_method2, n_simulations, rep(gamma1_Poisson,n_simulations)))
  sigma2iid_perf_indexes_PE_method2 = unlist(my_bias_RMSE_coverage(posterior_sigma2iid_PE_method2, n_simulations, rep(sigma2iid_Poisson,n_simulations)))
  
  
  #--- Export Poisson parameter performance
  
  output_Poissonpar_performance = rbind(
    c(gamma0_perf_indexes_PI_method1,
      gamma1_perf_indexes_PI_method1,
      sigma2iid_perf_indexes_PI_method1),
    c(gamma0_perf_indexes_PI_method2,
      gamma1_perf_indexes_PI_method2,
      sigma2iid_perf_indexes_PI_method2),
    c(gamma0_perf_indexes_FF_method1,
      gamma1_perf_indexes_FF_method1,
      sigma2iid_perf_indexes_FF_method1),
    c(gamma0_perf_indexes_FF_method2,
      gamma1_perf_indexes_FF_method2,
      sigma2iid_perf_indexes_FF_method2),
    c(gamma0_perf_indexes_PE_method1,
      gamma1_perf_indexes_PE_method1,
      sigma2iid_perf_indexes_PE_method1),
    c(gamma0_perf_indexes_PE_method2,
      gamma1_perf_indexes_PE_method2,
      sigma2iid_perf_indexes_PE_method2)
    )
  
  
  output_Poissonpar_performance = data.frame(method_exp = c(rep(1:2,3)),
                                             method_pois = c(rep(c("PI","FF","PE"),each=2)),
                                             scenario = rep(scenario_n,6),
                                             output_Poissonpar_performance)
  colnames(output_Poissonpar_performance)[-c(1:3)] = do.call(paste0,
                                                        expand.grid(c("bias.","rmse.","coverage."),c("gamma0","gamma1","sigma2iid")))

  
  
  if(!file.exists(paste0(getwd(),"/Output/Poisson_",name_file_performance))) {
    write.table(output_Poissonpar_performance,
                paste0(getwd(),"/Output/Poisson_",name_file_performance),
                sep=";",row.names = F)
  } else {
    write.table(output_Poissonpar_performance,
                paste0(getwd(),"/Output/Poisson_",name_file_performance),
                sep=";",row.names = F, col.names = F, append=T)
  }
  
  ##############################################################
  # Plot posterior distributions Poisson parameters
  ##############################################################
  #gamma1
  pdf(paste0(getwd(),"/Figures/","Posteriors_Poisson_gamma1_scenarioN",scenario_n,".pdf"))
  
  m.posterior_gamma1_PI_method1 = matrix(posterior_gamma1_PI_method1,ncol=n_simulations)
  m.posterior_gamma1_PI_method2 = matrix(posterior_gamma1_PI_method2,ncol=n_simulations)
  m.posterior_gamma1_FF_method1 = matrix(posterior_gamma1_FF_method1,ncol=n_simulations)
  m.posterior_gamma1_FF_method2 = matrix(posterior_gamma1_FF_method2,ncol=n_simulations)
  m.posterior_gamma1_PE_method1 = matrix(posterior_gamma1_PE_method1,ncol=n_simulations)
  m.posterior_gamma1_PE_method2 = matrix(posterior_gamma1_PE_method2,ncol=n_simulations)
  
  #par(mfrow=c(2,2))
  #Values ll together across iterations
  density_gamma1_PI_method1 = density(posterior_gamma1_PI_method1)
  density_gamma1_PI_method2 = density(posterior_gamma1_PI_method2)
  density_gamma1_FF_method1 = density(posterior_gamma1_FF_method1)
  density_gamma1_FF_method2 = density(posterior_gamma1_FF_method2)
  density_gamma1_PE_method1 = density(posterior_gamma1_PE_method1)
  density_gamma1_PE_method2 = density(posterior_gamma1_PE_method2)
  
  plot(density_gamma1_PI_method1,
       xlim=range(c(density_gamma1_PI_method1$x,density_gamma1_FF_method1$x,density_gamma1_PE_method1$x,density_gamma1_PE_method2$x)),
       ylim=range(c(density_gamma1_PI_method1$y,density_gamma1_PI_method2$y,
                    density_gamma1_FF_method1$y,density_gamma1_FF_method2$y,
                    density_gamma1_PE_method1$y,density_gamma1_PE_method2$y)),
       ylab = expression(gamma[1]),main="")
  lines(density_gamma1_PI_method2,col="red")
  
  lines(density_gamma1_FF_method1,lty=2,lwd=1)
  lines(density_gamma1_FF_method2,lty=2,lwd=1,col=2)
  
  lines(density_gamma1_PE_method1,lty=3,lwd=1)
  lines(density_gamma1_PE_method2,lty=3,lwd=1,col=2)
  
  abline(v=gamma1_Poisson,lty=4,col="green",lwd=2)
  #abline(v=median(posterior_gamma1_PI))
  #abline(v=median(posterior_gamma1_FF),lty=2,lwd=1)
  #abline(v=median(posterior_gamma1_PE),lty=3,lwd=1)
  legend("topleft",lty=c(1,1,2,2,3,3),
         legend=c("PI method1" , "PI method2", "FF method1", "FF method2", "PE method1", "PE method2"),
         col=c(1,2,1,2,1,2),box.lty=0)
  
  # # Single plot for each iteration
  # for(jj in 1:n_simulations){
  #   density_gamma1_PI = density(m.posterior_gamma1_PI[,jj])
  #   density_gamma1_FF = density(m.posterior_gamma1_FF[,jj])
  #   density_gamma1_PE = density(m.posterior_gamma1_PE[,jj])
  #   plot(density_gamma1_PI,
  #        xlim=range(c(density_gamma1_PI$x,density_gamma1_FF$x,density_gamma1_PE$x)),
  #        ylim=range(c(density_gamma1_PI$y,density_gamma1_FF$y,density_gamma1_PE$y)),
  #        ylab = expression(gamma[1]),main=jj)
  #   lines(density_gamma1_FF,lty=2,lwd=1)
  #   lines(density_gamma1_PE,lty=3,lwd=1)
  #   abline(v=gamma1_Poisson,lty=4,col=2,lwd=2)
  #   abline(v=median(m.posterior_gamma1_PI[,jj]))
  #   abline(v=median(m.posterior_gamma1_PE[,jj]),lty=3,lwd=1)
  #   legend("topleft",lty=c(1,2,3),legend=c("PI","FF","PE"),lwd=c(1,1,1),
  #          col=c(1,1,1),box.lty=0)
  # }
  dev.off()
  
  #gamma0
  pdf(paste0(getwd(),"/Figures/","Posteriors_Poisson_gamma0_scenarioN",scenario_n,".pdf"))
  
  m.posterior_gamma0_PI_method1 = matrix(posterior_gamma0_PI_method1,ncol=n_simulations)
  m.posterior_gamma0_PI_method2 = matrix(posterior_gamma0_PI_method2,ncol=n_simulations)
  m.posterior_gamma0_FF_method1 = matrix(posterior_gamma0_FF_method1,ncol=n_simulations)
  m.posterior_gamma0_FF_method2 = matrix(posterior_gamma0_FF_method2,ncol=n_simulations)
  m.posterior_gamma0_PE_method1 = matrix(posterior_gamma0_PE_method1,ncol=n_simulations)
  m.posterior_gamma0_PE_method2 = matrix(posterior_gamma0_PE_method2,ncol=n_simulations)
  
  #par(mfrow=c(2,2))
  #Values ll together across iterations
  density_gamma0_PI_method1 = density(posterior_gamma0_PI_method1)
  density_gamma0_PI_method2 = density(posterior_gamma0_PI_method2)
  density_gamma0_FF_method1 = density(posterior_gamma0_FF_method1)
  density_gamma0_FF_method2 = density(posterior_gamma0_FF_method2)
  density_gamma0_PE_method1 = density(posterior_gamma0_PE_method1)
  density_gamma0_PE_method2 = density(posterior_gamma0_PE_method2)
  
  plot(density_gamma0_PI_method1,
       xlim=range(c(density_gamma0_PI_method1$x,density_gamma0_FF_method1$x,density_gamma0_PE_method1$x,density_gamma0_PE_method2$x)),
       ylim=range(c(density_gamma0_PI_method1$y,density_gamma0_PI_method2$y,
                    density_gamma0_FF_method1$y,density_gamma0_FF_method2$y,
                    density_gamma0_PE_method1$y,density_gamma0_PE_method2$y)),
       ylab = expression(gamma[0]),main="")
  lines(density_gamma0_PI_method2,col="red")
  
  lines(density_gamma0_FF_method1,lty=2,lwd=1)
  lines(density_gamma0_FF_method2,lty=2,lwd=1,col=2)
  
  lines(density_gamma0_PE_method1,lty=3,lwd=1)
  lines(density_gamma0_PE_method2,lty=3,lwd=1,col=2)
  
  abline(v=gamma0_Poisson,lty=4,col="green",lwd=2)
  #abline(v=median(posterior_gamma0_PI))
  #abline(v=median(posterior_gamma0_FF),lty=2,lwd=1)
  #abline(v=median(posterior_gamma0_PE),lty=3,lwd=1)
  legend("topleft",lty=c(1,1,2,2,3,3),
         legend=c("PI method1" , "PI method2", "FF method1", "FF method2", "PE method1", "PE method2"),
         col=c(1,2,1,2,1,2),box.lty=0)
  
  dev.off()
  
  #sigma2iid
  pdf(paste0(getwd(),"/Figures/","Posteriors_Poisson_sigma2iid_scenarioN",scenario_n,".pdf"))
  
  m.posterior_sigma2iid_PI_method1 = matrix(posterior_sigma2iid_PI_method1,ncol=n_simulations)
  m.posterior_sigma2iid_PI_method2 = matrix(posterior_sigma2iid_PI_method2,ncol=n_simulations)
  m.posterior_sigma2iid_FF_method1 = matrix(posterior_sigma2iid_FF_method1,ncol=n_simulations)
  m.posterior_sigma2iid_FF_method2 = matrix(posterior_sigma2iid_FF_method2,ncol=n_simulations)
  m.posterior_sigma2iid_PE_method1 = matrix(posterior_sigma2iid_PE_method1,ncol=n_simulations)
  m.posterior_sigma2iid_PE_method2 = matrix(posterior_sigma2iid_PE_method2,ncol=n_simulations)
  
  #par(mfrow=c(2,2))
  #Values ll together across iterations
  density_sigma2iid_PI_method1 = density(posterior_sigma2iid_PI_method1)
  density_sigma2iid_PI_method2 = density(posterior_sigma2iid_PI_method2)
  density_sigma2iid_FF_method1 = density(posterior_sigma2iid_FF_method1)
  density_sigma2iid_FF_method2 = density(posterior_sigma2iid_FF_method2)
  density_sigma2iid_PE_method1 = density(posterior_sigma2iid_PE_method1)
  density_sigma2iid_PE_method2 = density(posterior_sigma2iid_PE_method2)
  
  plot(density_sigma2iid_PI_method1,
       xlim=range(c(density_sigma2iid_PI_method1$x,density_sigma2iid_FF_method1$x,density_sigma2iid_PE_method1$x,density_sigma2iid_PE_method2$x)),
       ylim=range(c(density_sigma2iid_PI_method1$y,density_sigma2iid_PI_method2$y,
                    density_sigma2iid_FF_method1$y,density_sigma2iid_FF_method2$y,
                    density_sigma2iid_PE_method1$y,density_sigma2iid_PE_method2$y)),
       ylab = expression(sigma^2),main="")
  lines(density_sigma2iid_PI_method2,col="red")
  
  lines(density_sigma2iid_FF_method1,lty=2,lwd=1)
  lines(density_sigma2iid_FF_method2,lty=2,lwd=1,col=2)
  
  lines(density_sigma2iid_PE_method1,lty=3,lwd=1)
  lines(density_sigma2iid_PE_method2,lty=3,lwd=1,col=2)
  
  abline(v=sigma2iid_Poisson,lty=4,col="green",lwd=2)
  #abline(v=median(posterior_sigma2iid_PI))
  #abline(v=median(posterior_sigma2iid_FF),lty=2,lwd=1)
  #abline(v=median(posterior_sigma2iid_PE),lty=3,lwd=1)
  legend("topleft",lty=c(1,1,2,2,3,3),
         legend=c("PI method1" , "PI method2", "FF method1", "FF method2", "PE method1", "PE method2"),
         col=c(1,2,1,2,1,2),box.lty=0)
  
  dev.off()  
  
  
  
}  

     