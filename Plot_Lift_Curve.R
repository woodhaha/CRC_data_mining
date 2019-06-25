LiftPlot = function(dataset) {
	lift <- data.frame(group = dataset$group)
	lift$logisreg <- predict(logisreg, dataset, type = "prob")[, "NonTumor"]
	lift$elastnet <- predict(elastnet, dataset, type = "prob")[, "NonTumor"]
	lift$RF <- predict(RF, dataset, type = "prob")[, "NonTumor"]
	lift$SVM <- predict(SVM, dataset, type = "prob")[, "NonTumor"]
	lift$NNET <- predict(NNET, dataset, type = "prob")[, "NonTumor"]
	lift_obj <- lift(group ~ logisreg + elastnet + RF + SVM + NNET, data = lift)
	ggplot(lift_obj, values = 60) + theme_bw() + theme(legend.position = "top")+
	labs(title = "Lift Curves of Models on dataset")+ 
	theme(plot.title = element_text(size = 15, hjust = 0.5, lineheight = 0.2))
}
