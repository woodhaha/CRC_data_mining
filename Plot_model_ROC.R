Plot_model_ROC <- function(models, dataset, color = "#1C61B6", add = F) {
	par(mar = c(4, 4, 6, 4))
	par(mfrow = c(1, 1))
	model.pred <- predict(models, newdata = dataset, type = "prob")
	model.roc <- roc(dataset$group, model.pred$Tumor)
	plot(model.roc, col = color, grid = c(0.1, 0.1), las = 1,
	   print.auc = F, print.auc.x = 1, print.auc.y = 0.4, print.auc.cex = 1, print.auc.col = "red",
		print.thres = T, print.thres.adj = c(-.05, 1.25), print.thres.cex = 0.8, add = add, main = "ROC Curve of Model")
 }