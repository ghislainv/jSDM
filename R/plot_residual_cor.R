## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

## Plot of residual correlation matrix (posterior mean estimator, and reordered by first principal component)
plot_residual_cor <- function(mod) {
	lv2.cor <- get_residual_cor(mod)
	lv2.cor$reorder.cor.mean <- corrplot::corrMatOrder(lv2.cor$cor.mean, order = "FPC", hclust.method = "average")
	rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- colnames(mod$model_spec$presences)
	
	par(cex=1, cex.main=1.5)
	corrplot::corrplot(lv2.cor$cor.mean[lv2.cor$reorder.cor.mean,lv2.cor$reorder.cor.mean], diag = F,
			 type = "lower", title = "Residual Correlation Matrix from LVM", mar = c(1,1,3,1),
			 method = "color", tl.srt = 45, tl.cex = 0.5)
}

# End