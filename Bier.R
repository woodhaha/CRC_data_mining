Bier = function(dataset) {
	require(prodlim)
	require(pec)
	require("doMC")
	require(riskRegression)
	registerDoMC()
	pecdata = dataset[, colnames(dataset) %in% c(cox_vip_list, "stage", "OS", "status")]
	Models <- list(
				"CoxPH.stage" = coxph(Surv(OS, status) ~ stage, data = pecdata, x = T, y = T),
				"CoxPH" = coxph(Surv(OS, status) ~ ., data = pecdata, x = T, y = T),
                "RandomForest" = rfsrc)
	Bier.imp <- pec::pec(Models, formula =Hist(OS, status) ~., data = pecdata,
	cens.model = "marginal", splitMethod = "bootcv", M = round(nrow(pecdata)*0.6), B =100, keep.index = T,
	multiSplitTest = T, confInt = T, confLevel = 0.95, exact = T, verbose = T, maxtime = 3000 , eval.times = (seq(0, 3650, 30)))
	return(Bier.imp)
}