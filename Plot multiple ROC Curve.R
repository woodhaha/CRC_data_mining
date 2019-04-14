Plot_multi_ROC <- function(x){
logisreg.pred <- predict(logisreg, newdata = x, type = "prob")
logisreg.roc <- roc(x$group, logisreg.pred$Tumor)

plot(logisreg.roc,col="#1C61B6",
     print.auc = F, print.auc.x = 1, las=1,print.auc.y = 0.4, print.auc.cex = 1, print.auc.col = "red",
     print.thres =F, print.thres.adj = c(-.05, 1.25), print.thres.cex = 0.8,grid=c(0.1, 0.1),xlim=c(1,0),
     main = "ROC Curve of different Models")


elastnet.pred <- predict(elastnet, newdata =x, type = "prob")
elastnet.roc <- roc(x$group, elastnet.pred$Tumor)
plot(elastnet.roc, col = "#FC4E07",
     print.auc = F, print.auc.x = 1, print.auc.y = 0.4, print.auc.cex = 1, print.auc.col = "red",
     print.thres =F, print.thres.adj = c(-.05, 1.25), print.thres.cex = 0.8, add = T)


RF.pred <- predict(RF, newdata = x, type = "prob")
RF.roc <- roc(x$group, RF.pred$Tumor)
plot(RF.roc, col = "#00AFBB",
     print.auc = F, print.auc.x = 1, print.auc.y = 0.4, print.auc.cex = 1, print.auc.col = "red",
     print.thres =F, print.thres.adj = c(-.05, 1.25), print.thres.cex = 0.8,
      add = T)

SVM.pred <- predict(SVM, newdata = x, type = "prob")
SVM.roc <- roc(x$group, SVM.pred$Tumor)
plot(SVM.roc, col = "#008600",
     print.auc = F, print.auc.x = 1, print.auc.y = 0.4, print.auc.cex = 1, print.auc.col = "red",
     print.thres =F, print.thres.adj = c(-.05, 1.25), print.thres.cex = 0.8,
     add = T)


NNET.pred <- predict(NNET, newdata = x, type = "prob")
NNET.roc <- roc(x$group, NNET.pred$Tumor)
plot(NNET.roc, col = "#E7B800",
     print.auc = F, print.auc.x = 1, print.auc.y = 0.4, print.auc.cex = 1, print.auc.col = "red",
     print.thres =F, print.thres.adj = c(-.05, 1.25), print.thres.cex = 0.8, add = T)

legend("bottomright", legend = c(paste("Logistic Regression AUC: ", round(logisreg.roc[["auc"]], 3)),
								paste("Elastic Net AUC: ", round(elastnet.roc[["auc"]], 3)),
					            paste("Random Forest AUC: ", round(RF.roc[["auc"]], 3)),
					            paste("Support Vector Machine AUC: ", round(SVM.roc[["auc"]], 3)),
								paste("Neural Network AUC: ", round(NNET.roc[["auc"]], 3)))
								,col = c("#1C61B6", "#FC4E07", "#00AFBB", "#008600", "#E7B800"), bty = "n",
									cex = 0.6, lwd = 2, title = "Model performace", title.adj = 0.5, adj = c(0, 0))
}



