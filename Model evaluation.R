## Model evaluation using Accuracy(ACC),Area under curve (AUC) 
model_pred <- function(models, dataset) {
	options(digits = 4)
	predsclas <- predict(models, newdata = dataset)
	res = confusionMatrix(data = predsclas, reference = dataset$group, positive = "Tumor", mode = "prec_recall")
	dat = data.frame(
		                 Accuracy = res[["overall"]][["Accuracy"]],
				 Lower = res[["overall"]][["AccuracyLower"]],
				 Upper = res[["overall"]][["AccuracyUpper"]],
				 AccPVal = res[["overall"]][["AccuracyPValue"]],
				 Kappa = res[["overall"]][["Kappa"]],
				 McnemarPVal = res[["overall"]][["McnemarPValue"]],
				 Precision=res[["byClass"]][["Precision"]],
		                 Recall=res[["byClass"]][["Recall"]],
		                 F1=res[["byClass"]][["F1"]])
	return(dat)
}

model_ROC <- function(models, dataset) {
	options(digits = 4)
	pred <- predict(models, newdata = dataset, type = "prob")
	roc <- roc(dataset$group, pred$Tumor)
	return(as.numeric(auc(roc)))
}


model_evaluation <- function(models) {
	options(digits = 4)
	temp = rbind(model_pred(models, train), model_pred(models, test))
	temp$DataSets = c("Training", "Testing")
	temp$AUC = c(model_ROC(models, train), model_ROC(models, test))
	temp = dplyr::select(temp, DataSets, AUC, Accuracy, Lower, Upper, AccPVal,Precision,Recall,F1, Kappa, McnemarPVal)
	return(temp)
}


