#' deg2rad
#' degree to rad
#' @param deg degree
#' @noRd
deg2rad = function(deg) {(deg * pi) / (180)}

#' rad2deg
#' rad to degree
#' @param rad rad
#' @noRd
rad2deg = function(rad) {(rad * 180) / (pi)}


ff = function(x){(x-min(x))/(max(x)-min(x))}

#' curve_text
#' plot curved text
#' @param pos position in degree
#' @param label label
#' @param radius radius
#' @param reverse in reverse order
#' @param middle text in the middle
#' @param extend extend char lengths, default 1.1
#' @param ... graphics::text
#' @noRd
curve_text = function(pos = 0, label = "", radius = 5.0, reverse = FALSE,middle = FALSE,extend = 1.1, ...){
  # inspired by plotrix package
  chars = strsplit(label, split = "")[[1]]
  char_lens = graphics::strwidth(chars)*extend
  char_angles = char_lens / radius
  changrang = range(char_angles)
  char_angles[char_angles < changrang[2]/2] = changrang[2]/2
  
  if(middle & reverse) pos = pos - rad2deg(sum(char_angles)/2)
  if(middle & !reverse) pos = pos + rad2deg(sum(char_angles)/2)
  
  if(reverse) {
    angles = c(deg2rad(pos), deg2rad(pos)+cumsum(char_angles)[-length(chars)])
    angles = angles + char_angles/2
  } else {
    angles = c(deg2rad(pos), deg2rad(pos)-cumsum(char_angles)[-length(chars)])
    angles = angles - char_angles/2
  }
  
  for(i in 1:length(chars)) graphics::text(label = chars[i],
                                           x = cos((angles[i]))*(radius),
                                           srt = rad2deg(angles[i]) - 90+ 180*reverse,
                                           y = sin((angles[i]))*(radius),
                                           xpd = NA, adj = c(0.5, 0.5),...)
  return(max(angles))
  
}


#' add_curve
#' curve plotting 'engine'
#' @param p1 first point
#' @param p2 second points
#' @param n number of points for spline
#' @param spar smoothing value
#' @param col curve's color
#' @param species draw species line
#' @param radius radius
#' @param lwd curve lwd
#' @noRd
add_curve = function(p1 = NULL, p2 = NULL, n = 10, spar = 0.7, col = "black", species = TRUE, radius = 5.0, lwd = 1.0) {
  xxs1 = cos(deg2rad(p1[3]))* seq(0, radius, length.out = n)
  xxs2 = cos(deg2rad(p2[3]))* seq(0, radius, length.out = n)
  yys1 = sin(deg2rad(p1[3]))* seq(0, radius, length.out = n)
  yys2 = sin(deg2rad(p2[3]))* seq(0, radius, length.out = n)
  x = c(rev(xxs1), xxs2[-1])
  y = c(rev(yys1), yys2[-1])
  m = (p1[2] - p2[2])/(p1[1] - p2[1])
  a = rad2deg(atan(m))
  a = -(a+180)
  alpha = deg2rad(a)
  alpha2 = deg2rad(-a)
  rot = matrix(c(cos((alpha)), -sin((alpha)), sin((alpha)), cos((alpha))),2,2)
  rot2 = matrix(c(cos((alpha2)), -sin((alpha2)), sin((alpha2)), cos((alpha2))),2,2)
  tt = cbind(x,y) %*% rot
  sp = stats::smooth.spline(tt[,1], tt[,2],spar = spar,df = 6, w = c(10.0, rep(0.1,nrow(tt)-2), 10.0))
  tt2 = cbind(sp$x, sp$y)
  b = tt2 %*% rot2
  graphics::lines(b[,1], b[,2], col = col, lwd = lwd)
  
  x1 = c(cos(deg2rad(p1[3]))*(radius+0.1), cos(deg2rad(p1[3]))*(radius+0.3))
  x2 = c(cos(deg2rad(p2[3]))*(radius+0.1), cos(deg2rad(p2[3]))*(radius+0.3))
  y1 = c(sin(deg2rad(p1[3]))* (radius+0.1), sin(deg2rad(p1[3]))* (radius+0.3))
  y2 = c(sin(deg2rad(p2[3]))* (radius+0.1), sin(deg2rad(p2[3]))* (radius+0.3))
  if(species){
    graphics::segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = "darkgrey")
    graphics::segments(x0 = x2[1], x1 = x2[2], y0 = y2[1], y1 = y2[2],  col = "darkgrey")
  }
}


#' add_legend
#' add legend to circular plot
#'
#' @param cols colors for gradients
#' @param range gradient range
#' @param radius radius
#' @param angles angles, start and end values in degree
#' @noRd

add_legend = function(cols = 1:11, range = c(-1,1), radius = 5.0, angles = c(110, 70)){
  angles = seq(angles[1], angles[2], length.out = length(cols)+1)
  for(i in 2:length(angles)){
    xx1 = (radius+0.4)*cos( seq(deg2rad(angles[i-1]),deg2rad(angles[i]) ,length.out=50) )
    xx2 = (radius+0.7)*cos( seq(deg2rad(angles[i-1]),deg2rad(angles[i]) ,length.out=50) )
    yy1 = (radius+0.4)*sin( seq(deg2rad(angles[i-1]),deg2rad(angles[i]) ,length.out=50)  )
    yy2 = (radius+0.7)*sin( seq(deg2rad(angles[i-1]),deg2rad(angles[i]) ,length.out=50)  )
    graphics::polygon(c(xx1, rev(xx2)), c(yy1, rev(yy2)),border = NA, col = cols[i-1], xpd = NA)
    if(i == 2 || i == length(angles)) {
      if(i ==2) label = range[1]
      else label = range[2]
      tmp_a = (angles[i-1]+angles[i])/2
      graphics::text(srt = tmp_a-90,
                     x = (radius+0.99)*cos(deg2rad(tmp_a)),
                     y =  (radius+0.99)*sin(deg2rad(tmp_a)),
                     xpd = NA, labels = label)
    }
  }
}

#' add_species_arrows
#' add species arrows to circle
#' @param radius radius
#' @param label label between arrows
#' @param reverse reverse label
#' @param start start point for arrow in degree
#' @param end end point for arrow in degree
#' @noRd

add_species_arrows = function(radius = 5.0, label = "Species", reverse = TRUE, start = 150, end = 270) {
  
  # first
  angles = seq(150,195,length.out = 100)
  xx = cos(deg2rad(angles))*(radius+0.6)
  yy = sin(deg2rad(angles))*(radius+0.6)
  graphics::lines(xx, yy, xpd = NA)
  end = curve_text(195, label,radius = radius*1.12,reverse = reverse)
  # second
  angles = seq(rad2deg(end)+3,rad2deg(end)+45+8,length.out = 100)
  xx = cos(deg2rad(angles))*(radius*1.12)
  yy = sin(deg2rad(angles))*(radius*1.12)
  graphics::lines(xx, yy, xpd = NA)
  arrow_angle = max(angles)-2.8
  graphics::polygon(x = c(cos(deg2rad(arrow_angle))*(radius*1.10), cos(deg2rad(arrow_angle))*(radius*1.14), cos(deg2rad(max(angles)))*(radius*1.12), cos(deg2rad(arrow_angle))*(radius*1.10)),
                    y = c(sin(deg2rad(arrow_angle))*(radius*1.10), sin(deg2rad(arrow_angle))*(radius*1.14), sin(deg2rad(max(angles)))*(radius*1.12), sin(deg2rad(arrow_angle))*(radius*1.10)),col = "black", xpd = NA)
}


#' plot_associations
#' plot species-species associations
#'
#' @param R  matrix of correlation \eqn{R}
#' @param radius circle's radius
#' @param main title
#' @param circleBreak circle break or not
#' @param top number of top negative and positive associations to consider
#' @param occ species occurence data
#' @param env_effect environmental species effects \eqn{\beta} 
#' @param cols_association color gradient for association lines
#' @param cols_occurrence color gradient for species 
#' @param cols_env_effect color gradient for environmental effect 
#' @param lwd_occurrence lwd for occurrence lines
#' @param species_order order species according to : \tabular{ll}{
#'  \code{"abundance"} \tab their mean abundance at sites by default) \cr
#'  \code{"frequency"} \tab the number of sites where they occur \cr
#'  \code{"main environmental effect"} \tab their most important environmental coefficients \cr}
#' @param species_indices indices for sorting species
#' @details  After fitting the jSDM with latent variables, the \bold{fullspecies residual correlation matrix} : \eqn{R=(R_ij) avec i=1,\ldots, nspecies et j=1,\ldots, nspecies}{R=(R_ij) avec i=1,..., nspecies et j=1,..., nspecies}
#'  can be derived from the covariance in the latent variables such as : 
#' \tabular{lll}{
#' \eqn{\Sigma_{ij}}{Sigma_ij} \tab \eqn{= \lambda_i .\lambda_j' + 1}{= \lambda_i . \lambda_j' + 1} \tab if i=j \cr
#'          \tab \eqn{= \lambda_i .\lambda_j'}{= \lambda_i . \lambda_j'} \tab else, \cr}
#' this function represents the correlations computed from covariances :
#'\deqn{R_{ij} = \frac{\Sigma_{ij}}{\sqrt{\Sigma_ii\Sigma _jj}}}{R_ij = Sigma_ij / sqrt(Sigma_ii.Sigma _jj)}.
#' @author \tabular{l}{
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>\cr
#' Jeanne Clément <jeanne.clement16@laposte.net>\cr }
#' @references \tabular{l}{
#' Pichler M. and Hartig F. (2020). “A new method for faster and more accurate inference of species associations from big community data”. \cr}
#' @seealso \code{\link{jSDM-package}} \code{\link{get_residual_cor}} \code{\link{jSDM_binomial_probit}} \code{\link{jSDM_binomial_probit_long_format}} \code{\link{jSDM_binomial_logit}}  \code{\link{jSDM_poisson_log}} 
#' @examples 
#' library(jSDM)
#' # frogs data
#' data(mites, package="jSDM")
#' # Arranging data
#' PA_mites <- mites[,1:35]
#' # Normalized continuous variables
#' Env_mites <- cbind(mites[,36:38], scale(mites[,39:40]))
#' colnames(Env_mites) <- colnames(mites[,36:40])
#' Env_mites <- as.data.frame(Env_mites)
#' # Parameter inference
#' # Increase the number of iterations to reach MCMC convergence
#' mod <- jSDM_poisson_log(# Response variable
#'                         count_data=PA_mites,
#'                         # Explanatory variables
#'                         site_formula = ~  water + topo + density,
#'                         site_data = Env_mites,
#'                         n_latent=2,
#'                         site_effect="random",
#'                         # Chains
#'                         burnin=100,
#'                         mcmc=100,
#'                         thin=1,
#'                         # Starting values
#'                         alpha_start=0,
#'                         beta_start=0,
#'                         lambda_start=0,
#'                         W_start=0,
#'                         V_alpha=1,
#'                         # Priors
#'                         shape=0.5, rate=0.0005,
#'                         mu_beta=0, V_beta=10,
#'                         mu_lambda=0, V_lambda=10,
#'                         # Various
#'                         seed=1234, verbose=1)
#' # Calcul of residual correlation between species
#' R <- get_residual_cor(mod)$cor.mean
#' plot_associations(R, circleBreak = TRUE, occ = PA_mites, species_order="abundance")
#' # Average of MCMC samples of species enrironmental effect beta except the intercept
#' env_effect <- t(sapply(mod$mcmc.sp,colMeans)[grep("beta_",colnames(mod$mcmc.sp[[1]]))[-1],])
#' colnames(env_effect) <-  gsub("beta_", "", colnames(env_effect))
#' plot_associations(R, env_effect = env_effect, species_order="main env_effect")
#' @importFrom graphics polygon text 
#' @export
plot_associations = function(R, radius = 5.0, main = NULL, 
                             circleBreak = FALSE, top = 10L, occ = NULL, env_effect=NULL, 
                             cols_association = c("#FF0000", "#BF003F", "#7F007F", "#3F00BF", "#0000FF"),
                             cols_occurrence = c( "#BEBEBE", "#8E8E8E", "#5F5F5F", "#2F2F2F", "#000000"),
                             cols_env_effect =c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                                "#66A61E", "#E6AB02", "#A6761D", "#666666"),
                             lwd_occurrence = 1.0,
                             species_order="abundance",
                             species_indices = NULL
){
  
  if(species_order!="main env_effect"){
    ##### circle #####
    
    lineSeq = 0.94*radius
    nseg = 100
    graphics::plot(NULL, NULL, xlim = c(-radius,radius), ylim =c(-radius,radius),pty="s", axes = F, xlab = "", ylab = "")
    if(!is.null(main)) graphics::text(x = 0, y = radius*1.35, pos = 3, xpd = NA, labels = main)
    xx = lineSeq*cos( seq(0,2*pi, length.out=nseg))
    yy = lineSeq*sin( seq(0,2*pi, length.out=nseg))
    
    graphics::polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
    
    #### curves ####
    n = ncol(R)
    
    if(!is.null(species_indices))
      R = R[species_indices, species_indices]
    else {
      ### occ ###
      if(!is.null(occ)) {
        species_indices = order(apply(occ, 2, sum))
        R = R[species_indices, species_indices]
      }
    }
    
    sigmas = R[upper.tri(R)]
    upper = order(sigmas, decreasing = TRUE)[1:top]
    lower = order(sigmas, decreasing = FALSE)[1:top]
    cuts = cut(sigmas, breaks = seq(-1,1,length.out = length(cols_association) + 1))
    to_plot = (1:length(sigmas) %in% upper) | (1:length(sigmas) %in% lower)
    levels(cuts) = cols_association
    cuts = as.character(cuts)
    
    angles = seq(0,355,length.out = n+1)[1:(n)]
    xx = cos(deg2rad(angles))*lineSeq
    yy = sin(deg2rad(angles))*lineSeq
    counter = 1
    coords = cbind(xx, yy, angles)
    for(i in 1:n) {
      for(j in i:n){
        if(i!=j) {
          #if(to_plot[counter]) add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, lineSeq = lineSeq)
          if(to_plot[counter]) add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, radius = lineSeq)
          counter = counter + 1
          #cat(counter, "\n")
        }
      }
    }
    
    
    ### occ ###
    if(!is.null(occ)) {
      lineSeq = radius
      occ_abs = sort(apply(occ, 2, sum))
      occ_logs = log(occ_abs)
      cuts = cut(occ_logs, breaks = length(cols_occurrence))
      cols = cols_occurrence #colfunc(5)
      levels(cuts) = cols_occurrence
      for(i in 1:length(occ_logs)){
        p1 = coords[i,]
        x1 = c(cos(deg2rad(p1[3]))*(lineSeq+0.1), cos(deg2rad(p1[3]))*(lineSeq+0.3))
        y1 = c(sin(deg2rad(p1[3]))* (lineSeq+0.1), sin(deg2rad(p1[3]))* (lineSeq+0.3))
        graphics::segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = as.character(cuts[i]), lend = 1, lwd = lwd_occurrence)
      }
      add_legend(cols_association, angles = c(140,110), radius = radius)
      graphics::text(cos(deg2rad(123))*(lineSeq+0.7), sin(deg2rad(123))*(lineSeq*1.14), labels = "correlation", pos = 2, xpd = NA)
      add_legend(cols = cols, range = c(min(occ_abs), max(occ_abs)), angles = c(70,40), radius = radius)
      graphics::text(cos(deg2rad(53))*(lineSeq+0.7), sin(deg2rad(55))*(lineSeq*1.14), labels =paste0("Sp.",species_order), pos = 4, xpd = NA)
    }
    
    ### arrows
    
    if(isTRUE(circleBreak)) {
      graphics::segments(x0 = cos(deg2rad(-1))*(lineSeq*0.96), x1 = cos(deg2rad(-1))*(lineSeq*1.18),
                         y0 = sin(deg2rad(-1))*(lineSeq*0.96), y1 = sin(deg2rad(-1))*(lineSeq*1.18), xpd = NA)
      graphics::segments(x0 = cos(deg2rad(356))*(lineSeq*0.96), x1 = cos(deg2rad(356))*(lineSeq*1.18),
                         y0 = sin(deg2rad(356))*(lineSeq*0.96), y1 = sin(deg2rad(356))*(lineSeq*1.18), xpd = NA)
    }
  }
  
  if(species_order=="main env_effect"){
    # Species effect beta 
    effects = env_effect
    np <- ncol(effects)
    within = apply(effects, 1, function(e) which.max(abs(e)))
    max_eff = sapply(1:length(within), function(i) effects[i,within[i]])
    effect_comb = cbind(max_eff, within)
    effect_comb_ind = order(within, max_eff)
    R = R[effect_comb_ind, effect_comb_ind]
    sigmas = R[upper.tri(R)]
    upper = order(sigmas, decreasing = TRUE)[1:top]
    lower = order(sigmas, decreasing = FALSE)[1:top]
    cuts = cut(sigmas, breaks = seq(-1,1,length.out = length(cols_association) + 1))
    to_plot = 1:length(sigmas) %in% upper | 1:length(sigmas) %in% lower
    levels(cuts) = cols_association
    cuts = as.character(cuts)
    n = ncol(R)
    ##### circle #####
    lineSeq = 0.94*radius
    nseg = 100
    graphics::plot(NULL, NULL, xlim = c(-radius,radius), ylim =c(-radius,radius),pty="s", axes = F, xlab = "", ylab = "")
    if(!is.null(main)) graphics::text(x = 0, y = radius*1.35, pos = 3, xpd = NA, labels = main)
    xx = lineSeq*cos( seq(0,2*pi, length.out=nseg))
    yy = lineSeq*sin( seq(0,2*pi, length.out=nseg))
    polygon(xx,yy, col= "white", border = "black", lty = 1, lwd = 1)
    ## curves
    angles = seq(0,360,length.out = n+1)[1:(n)]
    xx = cos(deg2rad(angles))*lineSeq
    yy = sin(deg2rad(angles))*lineSeq
    counter = 1
    coords = cbind(xx, yy, angles)
    for(i in 1:n) {
      for(j in i:n){
        if(i!=j) {
          if(to_plot[counter]) add_curve(coords[i,], coords[j,], col = cuts[counter], n = 5, species = TRUE, radius = lineSeq)
          counter = counter + 1
        }
      }
    }
    cols = cols_env_effect
    d = 1.0
    effect_comb2 = effect_comb
    effect_comb2[,1] = ff(effect_comb[,1])
    effect_comb2 = cbind(effect_comb2, effect_comb[,1])
    for(i in sort(unique(within))) {
      sub = coords[within[effect_comb_ind] == i,]
      sub_eff = effect_comb2[effect_comb_ind,][effect_comb2[effect_comb_ind,2] == i,]
      from = sub[1,3]
      to = sub[nrow(sub),3]
      
      x = c((radius+d*(sub_eff[,1]))*cos(deg2rad(sub[,3]) ), 
            rev((radius+d/2)*cos(deg2rad(sub[,3]))))
      
      y = c((radius+d*(sub_eff[,1]))*sin(deg2rad(sub[,3])),
            rev((radius+d/2)*sin(deg2rad(sub[,3]))))
      
      colnames(env_effect)
      angleName = (from+to)/2
      if(angleName > 180) reverse = TRUE
      else reverse = FALSE
      curve_text(angleName, label = colnames(env_effect)[i],reverse = reverse, radius=radius+d, middle = TRUE, extend = 1.17, col = cols[i])
      #y
      polygon(x, y, xpd = NA,col = cols[i])
      text(srt = 0, 
           x = (radius+0.1+d)*cos(deg2rad(sub[1,3]+4)), 
           y =  (radius+0.1+d)*sin(deg2rad(sub[1,3]+4)), 
           xpd = NA, labels = round(min(sub_eff[,3]), 2), col = cols[i], cex = 0.9)
      
      text(srt = 0, 
           x = (radius+0.1+d)*cos(deg2rad(sub[nrow(sub),3]-4)), 
           y =  (radius+0.1+d)*sin(deg2rad(sub[nrow(sub),3]-4)), 
           xpd = NA, labels = round(max(sub_eff[,3]), 2), col = cols[i], cex = 0.9)
    }
    lineSeq <- radius+d+0.1
    add_legend(cols_association, angles = c(140,110), radius=lineSeq)
    graphics::text(cos(deg2rad(125))*(lineSeq+0.7), sin(deg2rad(125))*(lineSeq*1.14), labels = "correlation", pos = 2, xpd = NA)
  }
  
  controlCircular = list()
  controlCircular$radius = radius
  controlCircular$n = n
  return(invisible(controlCircular))
}
