CaliPlot = function(dataset) { 
	lift <- data.frame(group = dataset$group)
	lift$logisreg <- predict(logisreg, dataset, type = "prob")[, "NonTumor"]
        lift$elastnet <- predict(elastnet, dataset, type = "prob")[,"NonTumor"]
        lift$RF <- predict(RF, dataset, type = "prob")[,"NonTumor"]
        lift$SVM <- predict(SVM, dataset, type = "prob")[,"NonTumor"]
        lift$NNET <- predict(NNET, dataset, type = "prob")[,"NonTumor"]
      
	cal_obj <- calibration(group ~ logisreg + elastnet + RF+SVM+NNET, data = lift,cuts = 10)
        
	ggplot(cal_obj, type = "l", auto.key = list(columns = 3, lines = T, points = F)) +
	theme_bw() + theme(legend.position = "top") +
	labs(title = "Calibration Curves of Models on dataset") +
	theme(plot.title = element_text(size = 15, hjust = 0.5, lineheight = 0.2))
}
