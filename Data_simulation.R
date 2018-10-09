############################################################
# This function simulates the exposure data at the point level,
# computes the exposure area mean and simulates the Poisson data
# Last update: 01/02/2017
  ############################################################

Data_simulation = function(scenario_n,
                           sample.by,
                           prop.points.for.area,
                           res_pred_grid,
                           range.scale, sigma2_xi, b0, b1, sigma2_e,
                           z_sim_grid,
                           gamma1_Poisson, gamma0_Poisson,
                           sigma2iid_Poisson, fixed_E) {

  save.plot = T
  BW = F #black and white plots
  
  #--- Import the shapefile
  shapefile_sim = readShapePoly(system.file("etc/shapes/bhicv.shp",package="spdep")[1]) #coordinates in km 
  colnames(shapefile_sim@data)[2] = "AREAKEY" 
  
  #Transform the shapefile from km to m
  #library(rgdal)
  #shapefile_sim = spTransform(shapefile_sim, CRS("+proj=utm +ellps=WGS84 +datum=WGS84 +units=km +no_defs")) #proj=utm
  #writePolyShape(shapefile_sim,"ShapeFile_CCG/CCG_BSC Feb2013  (clipcoast 200m).shp")
  
  n_areas = nrow(shapefile_sim@data)
  shapefile_sim_centroids = gCentroid(shapefile_sim, byid = T)
  centroids_df = data.frame(AREAKEY = shapefile_sim@data$AREAKEY, centroids = shapefile_sim_centroids)
  colnames(centroids_df)[2:3] = c("x1", "x2")
  
  # require(rgeos)
  # plot(shapefile_sim)
  # centroids <- gCentroid(shapefile_sim, byid=TRUE)
  # centroidLons <- coordinates(centroids)[,1]
  # centroidLats <- coordinates(centroids)[,2]
  # text(centroidLons, centroidLats, labels=shapefile_sim@data$AREAKEY,  cex=.7)
  # text(centroidLons, centroidLats, labels=round(shapefile_sim@data$polygonAREA/sum(shapefile_sim@data$polygonAREA)*100,2),  cex=.7)
  
  ############################################################
  # Create a high resolution grid of points  for exposure simulation
  # The resolution is highres_simulation_reggrid
  ############################################################
  sim_grid_x1 = seq(shapefile_sim@bbox[1, 1] - 0.1, shapefile_sim@bbox[1, 2] + 0.1, l = highres_simulation_reggrid) 
  sim_grid_x2 = seq(shapefile_sim@bbox[2, 1] - 0.1, shapefile_sim@bbox[2, 2] + 0.1, l = highres_simulation_reggrid) 
  
  sim_grid = expand.grid(x = sim_grid_x1, y = sim_grid_x2)
  coordinates(sim_grid) <- ~ x + y
  colnames(sim_grid@coords) = c("x1", "x2")
  proj4string(sim_grid) = proj4string(shapefile_sim)
  
  # plot(shapefile_sim)
  # plot(sim_grid, add=T)
  
  ############################################################
  # Exposure simulation
  ############################################################
  
  #--- Associate each point to an area of the shapefile and select points inside the shapefile
  sim_grid_df = over(sim_grid, shapefile_sim)
  sim_grid_df$z = z_sim_grid
  sim_grid_df = cbind(sim_grid_df,sim_grid@coords)
  
  #--- Set the parameters for the exposure simulation
  max_dist = max(dist(as.matrix(sim_grid@coords)))
  #max(dist(as.matrix(shapefile_sim_centroids@coords)))
  range_xi = max_dist * range.scale
  kappa_xi = sqrt(8) / range_xi
  
  #--- Simulate the spatial process and then the exposure field
  #cat("--- I'm simulating the data ... it can take a while" , "\n")
  xi = simulate_GF(coords = sim_grid@coords,
                   kappa = kappa_xi,
                   variance = sigma2_xi)
  exposure = b0 + b1 * sim_grid_df$z + xi #+ rnorm(length(sim_grid),mean = 0,sd = sqrt(sigma2_e)) #with measurement error
  sim_grid_df$exposure = exposure
  
  ##############################################################
  # ------ Plot the values at the sim_grid sites
  ##############################################################
  #--- Just select the points inside the shapefile
  which_inside = which(!is.na(sim_grid_df$AREAKEY))
  sim_grid_inside_df = sim_grid_df[which_inside, ]
  sim_grid_inside = sim_grid[which_inside]
  #areas_with_nopoints = which(table(sim_grid_inside_df$AREAKEY)==0) #very small areas without points inside (to be exluded when sampling monitoring stations)
  
  cutoffs <- quantile(sim_grid_inside_df$exposure, seq(from = 0, to = 1, by = 0.25))
  exposure_classes <- cut(sim_grid_inside_df$exposure, breaks = cutoffs, include.lowest = TRUE)
  if (save.plot) {
    if (BW) {
      exposure_colors = colorRampPalette(c("white", "black"))(length(cutoffs))[as.numeric(exposure_classes)]
      legend_col = colorRampPalette(c("white", "black"))(length(cutoffs))
    } else{
      #x_colors = topo.colors(length(cutoffs))[as.numeric(exposure_classes)]
      x_colors = colorRampPalette(c("green", "red"))(4) [as.numeric(exposure_classes)]
      #legend_col = topo.colors(length(cutoffs))
      legend_col = colorRampPalette(c("green", "red"))(4) 
    }
    
    pdf(paste0(getwd(), "/Data/True_exposure.pdf"), width=5, height=6)
    # True exposure values at the point level
    par(mar = c())
    plot(shapefile_sim, las = 1, lwd = 2)
    points(sim_grid_inside, bg = x_colors, pch = 21, cex=0.5)
    legend(
      "topleft",
      pch = rep(21, length(cutoffs)),
      pt.cex = rep(2.5, length(cutoffs)),
      pt.bg = rev(legend_col),
      legend = rev(levels(exposure_classes)),
      bty = "n"
    )
    plot(shapefile_sim, add = T, lwd = 1)
    dev.off()
  }
  
  ############################################################
  # Compute the TRUE *AREA* exposure as mean of the sim_grid points inside each area
  ############################################################
  true_area_exposure = aggregate(exposure ~ AREAKEY,
                                 data=sim_grid_inside_df,
                                 FUN=mean)
  colnames(true_area_exposure)[2] = "true_area_exposure"

  #--- Divide the area values in classes (NB: the cutoffs are the same computed at the point level)
  true_area_exposure$true_area_exposure_class <- cut(true_area_exposure$true_area_exposure,
                                                breaks=cutoffs,include.lowest=TRUE)
  
  #--- Include the "true" area (mean) value in the shapefile
  shapefile_sim = merge(shapefile_sim,
                              true_area_exposure,
                              by = "AREAKEY",
                              all.x = T)
  if(save.plot){
    pdf(paste0(getwd(), "/Data/True_exposure_areas.pdf"), width=5, height=6)
    # True exposure values at the point level
    par(mar = c())
    p = spplot(shapefile_sim,
         "true_area_exposure_class",
         col.regions=legend_col,
         par.settings = list(axis.line = list(col = 'transparent')))
    print(p)
    dev.off()
  }
  
  #--- Compute the variability of the true area points inside each area
  true_area_exposure_var = aggregate(exposure ~ AREAKEY,
                                 data=sim_grid_inside_df,
                                 FUN=var)
  colnames(true_area_exposure_var)[2] = "var_true_area_exposure"
  shapefile_sim = merge(shapefile_sim,
                        true_area_exposure_var,
                        by = "AREAKEY",
                        all.x = T)
  
  ############################################################
  # Create a regular grid for the prediction
  # grid resolution = highres_simulation_reggrid/res_pred_grid
  ############################################################
  pred_grid_x1 = unique(sim_grid@coords[, 1])[seq(1, highres_simulation_reggrid,
                                                  by = highres_simulation_reggrid / res_pred_grid)]
  pred_grid_x2 = unique(sim_grid@coords[, 2])[seq(1, highres_simulation_reggrid,
                                                  by = highres_simulation_reggrid / res_pred_grid)]
  
  pred_grid = expand.grid(x = pred_grid_x1, y = pred_grid_x2)
  coordinates(pred_grid) <- ~ x + y
  proj4string(pred_grid) = proj4string(sim_grid)
  # Transform the regular grid into a SpatialPolygons object
  pred_grid_poly = as.SpatialPolygons.GridTopology(points2grid(pred_grid))
  proj4string(pred_grid_poly) = proj4string(shapefile_sim)
  #plot(shapefile_sim)
  #plot(pred_grid_poly,add=T,border="grey")
  
  # Change the name of the pixels of the prediction grid (from left to right from bottom to top)
  ncol_pred_grid = res_pred_grid
  nrow_pred_grid = res_pred_grid
  mydesiredseq = c()
  for (i in (nrow_pred_grid - 1):0) {
    temp = seq((ncol_pred_grid * i) + 1, l = ncol_pred_grid)
    mydesiredseq = c(mydesiredseq, temp)
  }
  for (i in 1:length(pred_grid_poly@polygons)) {
    pred_grid_poly@polygons[[i]]@ID = paste("g", mydesiredseq[i], sep = "")
  }
  
  if(save.plot){
  pdf(paste0(getwd(), "/Data/Regular_grid.pdf"),width=5, height=6)  
  plot(pred_grid_poly,border = 'blue')
  plot(shapefile_sim,add=T)
  points(pred_grid,col=2,pch = 21, cex=0.5)
  #text(coordinates(pred_grid_poly),sapply(pred_grid_poly@polygons, function(x) x@ID),cex=0.75)
  dev.off()
  }
  
  
  which_pred_grid_in_sim_grid = over(pred_grid,sim_grid)
  exposure_classes_pred <- cut( sim_grid_df$exposure[which_pred_grid_in_sim_grid], breaks = cutoffs, include.lowest = TRUE)
  x_colors_pred = colorRampPalette(c("green", "red"))(4) [as.numeric(exposure_classes_pred)]
  
  if(save.plot){
    pdf(paste0(getwd(), "/Data/True_exposure_predpoints.pdf"), width=5, height=6)
    
    plot(pred_grid,cex=0)
    #plot(shapefile_sim,add=T)
    #points(pred_grid$x[1:40],pred_grid$y[1:40],col=2)
    #points(sim_grid_df$x1[1:130],sim_grid_df$x2[1:130],col=3) 
    points(pred_grid, bg = x_colors_pred, pch = 21, cex=0.5)
    legend(
      "topleft",
      pch = rep(21, length(cutoffs)),
      pt.cex = rep(2.5, length(cutoffs)),
      pt.bg = rev(legend_col),
      legend = rev(levels(exposure_classes_pred)),
      bty = "n"
    )
    plot(shapefile_sim, add = T, lwd = 1)
    dev.off()
  }
  
 
  
  ########################################################################################################################
  #--- Compute the INTERSECTIONS between the shapefile and the grid (only if it doesn't already exist)
  ########################################################################################################################
  if (!file.exists(
    paste0(getwd(),"/Data","/weights_grid_areas_intersection_res_pred_grid_",res_pred_grid,".Rdata"))) {
    
    weights_grid_areas_intersection = data.frame()
    cat("--- I'm computing the intersection weight matrix ... it takes a while" , "\n")
    
    for (i in 1:length(shapefile_sim@data$AREAKEY)) {
      # cat("--- Area n.",i,"---- of", length(shapefile_sim@data$AREAKEY), "\n" )
      # Select an area
      sel_area = shapefile_sim[i, ]
      if (!gIsValid(sel_area)) {
        sel_area = gBuffer(sel_area, width = 0, byid = TRUE)
      }
      
      # Plot the selected area
      # plot(shapefile_sim)
      # plot(sel_area,col="orange", add=T,border="red", lwd=4)
      
      # Compute the intersection between the selected area and the grid
      inter = gIntersection(pred_grid_poly, sel_area, byid = TRUE)
      if (is.null(inter))
        next #Skip if there is no intersection between the shapefile and the grid (eg. Isles of Scilly)
      
      # Define a data frame with the name+area intersection pixels and the weights given as area/total_area
      inter_w = data.frame(
        pixel = sapply(inter@polygons, function(x) x@ID),
        area = sapply(inter@polygons, function(x) x@area))
      inter_w$weights = inter_w$area / sum(inter_w$area)
      inter_w$AREAKEY = sel_area@data$AREAKEY
      # Check the areas (total area of the polygon - area given by the inserctions --> should be close to 0)
      abs_diff_area = abs(sum(sapply(inter@polygons, function(x) x@area)) - sel_area@polygons[[1]]@area)
      
      if (abs_diff_area  > 0.0005) {
        warning(
          paste("----------- ATTENTION: there's something wrong with the areas! for area i=",i,"The abs_diff_area is ", round(abs_diff_area, 6)))
      }
      
      # Plot the intersection
      #plot(sel_area)
      #plot(pred_grid_poly,add=T,border = 'blue')
      
      #Just for the paper area 1
      #plot(pred_grid_poly,border = 'blue',
      #ylim=sel_area@bbox[2,]+c(0.3,-0.25),xlim=c(-44.4,-43.84))
      #plot(sel_area,add=T)
      
      #plot(inter, add=TRUE, density = 5, col = 'blue')
      #class(inter)
      #require(rgeos)
      #centroids <- gCentroid(inter, byid=TRUE)
      #centroidLons <- coordinates(centroids)[,1]
      #centroidLats <- coordinates(centroids)[,2]
      #points(centroidLons, centroidLats)
      #points(pred_grid,col=2,pch=19)
      #text(centroidLons, centroidLats, labels=round(inter_w$weights,3),  cex=1.2)
      
      #just for the paper
      #points(pred_grid@coords[substring(inter_w$pixel,2,4),],bg=2,pch=21)
      #text(pred_grid@coords[substring(inter_w$pixel,2,4),1],
      #pred_grid@coords[substring(inter_w$pixel,2,4),2]-0.025,
      #labels=round(inter_w$weights,3),  cex=1.2)
      
      weights_grid_areas_intersection = rbind(weights_grid_areas_intersection, inter_w)
    }
    
    save(file = paste0(getwd(),
                       "/Data",
                       "/weights_grid_areas_intersection_res_pred_grid_",res_pred_grid,".Rdata"),
         weights_grid_areas_intersection)
    
  } else {
    
    load(
      paste0(getwd(),"/Data",
             "/weights_grid_areas_intersection_res_pred_grid_",
             res_pred_grid,".Rdata"))
  }
  
  ########################################################################################################################
  #--- Compute the pred points inside each area or the closest one  (only if it doesn't already exist)
  ########################################################################################################################
  if (!file.exists(
    paste0(getwd(),"/Data","/weights_points_areas_inside_closest_",res_pred_grid,".Rdata"))) {
    
    weights_points_areas_inside_closest = data.frame()
    cat("--- I'm computing the inside points weight matrix ... it takes a while" , "\n")
    
    for (i in 1:length(shapefile_sim@data$AREAKEY)) {
      
      sel_area = shapefile_sim[i, ]
      #Pred grid points inside the area
      pointsinside = which(!is.na(over(pred_grid,sel_area)$AREAKEY))
      
      if(length(pointsinside)>0){
        points_df = data.frame(pixel = pointsinside,
                               weights = rep(1/length(pointsinside),length(pointsinside)),
                               AREAKEY = shapefile_sim@data$AREAKEY[i])
      } else {
        nearestpoints = apply(gDistance(pred_grid, gCentroid(sel_area), byid=TRUE), 1, which.min)
        points_df = data.frame(pixel = nearestpoints,
                               weights = rep(1/length(nearestpoints),length(nearestpoints)),
                               AREAKEY = shapefile_sim@data$AREAKEY[i])
        
      }
    
      #Just for the paper area 1
      #plot(pred_grid_poly,border = 'blue',
      #ylim=sel_area@bbox[2,]+c(0.3,-0.25),xlim=c(-44.4,-43.84))
      #plot(sel_area,add=T)
      #plot(inter, add=TRUE, density = 5, col = 'blue')
      #points(pred_grid@coords[points_df$pixel,],bg=2,pch=21)
      #text(pred_grid@coords[points_df$pixel,1],
      #pred_grid@coords[points_df$pixel,2]-0.025,
      #labels=round(points_df$weights,3),  cex=1.2)
      
      
      
        weights_points_areas_inside_closest = rbind(weights_points_areas_inside_closest, points_df)
    }
    
    save(file = paste0(getwd(),
                       "/Data",
                       "/weights_points_areas_inside_closest_",res_pred_grid,".Rdata"),
         weights_points_areas_inside_closest)
    
  } else {
    
    load(
      paste0(getwd(),"/Data",
             "/weights_points_areas_inside_closest_",
             res_pred_grid,".Rdata"))
  }
  
  
  ########################################################################################################################
  #--- Prediction dataset
  ########################################################################################################################
  which_in_pred_grid = over(pred_grid, sim_grid)
  sim_grid_inside_df_pred = sim_grid_df[which_in_pred_grid, ]
  
  
  #############################################################
  # Sample the monitoring station sites (they are fixed for each iteration of each scenario)
  ############################################################
  sim_grid_inside_df2 = data.frame(sim_grid_inside_df, sim_grid_inside@coords)
  sim_grid_inside_df2 = merge(sim_grid_inside_df2, true_area_exposure, by="AREAKEY")
  
  if (!file.exists(paste0(getwd(),
                          "/Data",
                          "/monitoringstationsdata_scenario",scenario_n,".Rdata"))) {

        
    #--- Select the areas according to the sampling method ("all areas")
    selected_area = shapefile_sim@data$AREAKEY 
      
    #Compute how many monitoring stations we have for each area (according to prop.points.for.area)
    npoints.for.area = round(table(sim_grid_inside_df2$AREAKEY) * prop.points.for.area,0)
    
    #--- Select randomly a given number of locations within each selected area
    sim_grid_inside_df2_reduced = sim_grid_inside_df2[which(sim_grid_inside_df2$AREAKEY %in% selected_area), ]
    sim_grid_inside_df2_splitted = split(sim_grid_inside_df2_reduced,
                                         list(sim_grid_inside_df2_reduced$AREAKEY))
    #just keep elements with nrow>0 (remove areas with no points inside)
    area.nopoints = which(table(sim_grid_df$AREAKEY) == 0) #area without grid points inside
    sim_grid_inside_df2_splitted = sim_grid_inside_df2_splitted[sapply(sim_grid_inside_df2_splitted, nrow) >0]
    
    #all the areas, sample a given % of points
    if (length(area.nopoints) > 0) { #if there are some areas without points inside
        samples <- lapply(1:length(sim_grid_inside_df2_splitted),
                          function(x) {
                            sim_grid_inside_df2_splitted[[x]][sample(1:nrow(sim_grid_inside_df2_splitted[[x]]),
                                                                     npoints.for.area[-area.nopoints][x],
                                                                     replace = FALSE), ]})
    } else {
        samples <- lapply(1:length(sim_grid_inside_df2_splitted),
                          function(x) {
                            sim_grid_inside_df2_splitted[[x]][sample(1:nrow(sim_grid_inside_df2_splitted[[x]]),
                                                                     npoints.for.area[x],
                                                                     replace = FALSE), ]})
    }
    
    
    #--- Define the df with the monitoring stations
    selected_loc_df <- do.call(rbind, samples)
    monitoringstations.coords = selected_loc_df[, c("x1", "x2")]
    monitoringstations.data = selected_loc_df
    #dim(monitoringstations.coords)
    #freqmonitoringstats = as.numeric(table(monitoringstations.data$AREAKEY))
    #summary(freqmonitoringstats)
    
    
        
    save(file = paste0(getwd(),
                       "/Data",
                       "/monitoringstationsdata_scenario",scenario_n,".Rdata"),
                monitoringstations.coords, monitoringstations.data)
    
    
    
    } else { #retrieve just the stations coordinates and update the exposure values
      
      load(paste0(getwd(),
                "/Data",
                "/monitoringstationsdata_scenario",scenario_n,".Rdata"))
    
      #-- Merge the dataframe with all the exposure simulation points (sim_grid_inside_df2)
      #-- with the selected monitoring stats dataframe containing just the x1, x2 coordinates
      #-- This is required to select again the same monitoring stations (obviously withe different exposure values)
      monitoringstations.data = merge(sim_grid_inside_df2,
                                       monitoringstations.data[, c("x1","x2")], #just x1, x2
                                       by=c("x1","x2"), sort=F)
      
      monitoringstations.coords = monitoringstations.data[,c("x1","x2")]
      
    }
  
  
  if (save.plot) {
    xstats_classes <-cut(monitoringstations.data$exposure,breaks = cutoffs,include.lowest = TRUE)
    
    if (BW) {
      x_colors = colorRampPalette(c("white", "black"))(length(cutoffs))[as.numeric(x_classes)]
      legend_col = colorRampPalette(c("white", "black"))(length(cutoffs))
    } else{
      #x_colors = topo.colors(length(cutoffs))[as.numeric(xstats_classes)]
      #legend_col = topo.colors(length(cutoffs))
      x_colors = colorRampPalette(c("green", "red"))(4) [as.numeric(xstats_classes)]
      legend_col = colorRampPalette(c("green", "red"))(4) 
    }
    
    pdf(paste0(getwd(),"/Data/Monitoringstats.pdf"), width=5, height=6)
    # True exposure values (x) at the point level
    par(mar = c())
    plot(shapefile_sim, las = 1)
    points(monitoringstations.coords, bg = x_colors, pch = 21, cex=0.5)
    legend(
      "topleft",
      pch = rep(21, length(cutoffs)),
      pt.cex = rep(2.5, length(cutoffs)),
      pt.bg = rev(legend_col),
      legend = rev(levels(xstats_classes)),
      bty = "n"
    )
    plot(shapefile_sim, add = T)
    dev.off()
  }
  
  ############################################################
  # Simulate the exposure for the monitoring stations (just add the measurement error)
  # monitoringstations.data$exposure does not change across iterations of a single scenario
  ############################################################
  meas_error = rnorm(length(monitoringstations.data$exposure),mean = 0,sd = sqrt(sigma2_e))
  monitoringstations.data$exposure_with_meas_err = monitoringstations.data$exposure + meas_error
    
  ############################################################
  # Simulate the Poisson data
  ############################################################
  log_rate_mean_exposure = gamma0_Poisson + gamma1_Poisson * shapefile_sim@data$true_area_exposure + rnorm(n_areas, sd = sqrt(sigma2iid_Poisson))
  y_Poisson_mean_exposure = rpois(n_areas, lambda = exp(log_rate_mean_exposure) * fixed_E)
  data_Poisson = data.frame("AREAKEY" = shapefile_sim@data$AREAKEY, "y_Poisson" = y_Poisson_mean_exposure)
  
  #Add Poisson simulated DATA to the shapefile
  shapefile_sim = merge(shapefile_sim, data_Poisson, by = "AREAKEY")
  
  #Set the expected values for all the areas and include them in the shapefile
  expected = rep(fixed_E, n_areas)
  data.expected = data.frame(expected = expected, "AREAKEY" = shapefile_sim@data$AREAKEY)
  shapefile_sim = merge(shapefile_sim, data.expected, by = "AREAKEY")
  
  #pdf(paste0(getwd(),"/Figures/Poisson/iteration_n",sim_i,"/y_expected_Poisson.pdf"))
  #par(mar=c())
  #print(spplot(shapefile_sim,c("expected","y_Poisson"),col.regions=topo.colors(16)))
  #dev.off()

  ############################################################
  # Save all the relevant data for next analysis
  ############################################################
  shapefile_sim_complete = shapefile_sim
  
  return(list(shapefile_sim_complete=shapefile_sim_complete,
              monitoringstations.data=monitoringstations.data,
              monitoringstations.coords=monitoringstations.coords,
              max_dist=max_dist,
              pred_grid_x1 = pred_grid_x1,
              pred_grid_x2 = pred_grid_x2,
              pred_grid=pred_grid,
              n.monitoringstats = nrow(monitoringstations.data),
              sim_grid_inside_df=sim_grid_inside_df,  
              sim_grid_inside_df_pred=sim_grid_inside_df_pred
              ))
  
  
###############################################################  
#stop!
###############################################################
  
  
  if(F){
  ###############################################################  
  ###############################################################
  # Compute intersection weights (between shapefile areas and regular grid)
  # this is required later for the INLA lincombs
  #
  # ATTENTION: this code can't run on the server due to some problems with rgeos package
  # Create locally the weights_grid_areas_intersection objection and then copy it on the server
  ###############################################################
  # Define for each area of the shapefile the intersection with the regular grid and compute the weights
  # This has to be repeated if you change the grid resolution 
  
  if(!cluster){
    weights_grid_areas_intersection = data.frame()
    for(i in 1:length(shapefile_sim@data$AREAKEY)){ 
      cat("--- Area n.",i,"---- of", length(shapefile_sim@data$AREAKEY), "\n" )
      # Select an area
      sel_area = shapefile_sim[i,]
      if(!gIsValid(sel_area)) {sel_area = gBuffer(sel_area,width=0,byid=TRUE)}
      
      # Plot the selected area
      #plot(shapefile_sim)
      #plot(sel_area,col="orange", add=T,border="red", lwd=4) 
      
      # Compute the intersection between the selected area and the grid
      inter = gIntersection(pred_grid_poly, sel_area, byid = TRUE)
      
      if(is.null(inter)) next #Skip if there is no intersection between the shapefile and the grid (eg. Isles of Scilly)
      
      # Define a data frame with the name+area intersection pixels and the weights given as area/total_area
      inter_w = data.frame(pixel = sapply(inter@polygons, function(x) x@ID),
                           area = sapply(inter@polygons, function(x) x@area))
      inter_w$weights = inter_w$area/sum(inter_w$area)
      inter_w$AREAKEY = sel_area@data$AREAKEY
      #inter_w$CODE = sel_area@data$CODE
      #inter_w$NAME = sel_area@data$NAME
      
      # Check the areas (total area of the polygon - area given by the inserctions --> should be close to 0)
      abs_diff_area = abs(sum(sapply(inter@polygons, function(x) x@area)) - sel_area@polygons[[1]]@area)
      if( abs_diff_area  > 0.0005){
        warning(paste("----------- ATTENTION: there's something wrong with the areas! for area i=",i,
                      "The abs_diff_area is ", round(abs_diff_area, 6)))
      }
      
      # Plot the intersection
      #plot(sel_area)
      #plot(pred_grid_poly,add=T,border = 'blue')
      #plot(inter, add=TRUE, density = 5, col = 'blue')
      #class(inter)
      #require(rgeos)
      #centroids <- gCentroid(inter, byid=TRUE)
      #centroidLons <- coordinates(centroids)[,1]
      #centroidLats <- coordinates(centroids)[,2]
      #text(centroidLons, centroidLats, labels=round(inter_w$weights,2),  cex=1.7)
      #text(centroidLons[1:7], centroidLats[1:7], labels=round(inter_w$weights[1:7],2),  cex=1.7)
      
      weights_grid_areas_intersection = rbind(weights_grid_areas_intersection, inter_w)
    }
    save(file=paste0(getwd(),"/Data","/weights_grid_areas_intersection_res_pred_grid_",res_pred_grid,".Rdata"),
         weights_grid_areas_intersection)  
  } else {
    load(paste0(getwd(),"/Data","/weights_grid_areas_intersection_res_pred_grid_",res_pred_grid,".Rdata")) 
  }
  } #close if(F)
}  
  

