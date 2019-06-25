Cindex <- function(dataset) {
	require(prodlim)
	require(pec)
	require("doMC")
	require(riskRegression)
	registerDoMC()
	cdata = dataset[, colnames(dataset) %in% c(cox_vip_list, "stage", "OS", "status")]
	Models <- list(
		      "CoxPH.stage" = coxph(Surv(OS, status) ~stage, data = cdata, x = T, y = T) ,
		      "CoxPH" = coxph(Surv(OS, status) ~., data =cdata, x = T, y = T),
	             "RandomForest" = rfsrc)
	Cimp=pec::cindex(Models, formula = Hist(OS, status)~., data = cdata,eval.times = (seq(0, 3650, 30)),
	splitMethod = "bootcv", M = round(nrow(cdata)*0.6), cens.model = "marginal",B=100,maxtime = 3000 )
	return(Cimp)
}
## cens.modelMethod for estimating inverse probability of censoring weigths:cox: A semi-parametric Cox proportional hazard model is fitted to the censoringtimesmarginal: The Kaplan-Meier estimator for the censoring timesnonpar: Nonparametric extension of the Kaplan-Meier for the censoring timesusing  symmetric  nearest  neighborhoods  –  available  for  arbitrary  many  stratavariables on the right hand side of argumentformulabut at most one continuousvariable.  See the documentation of the functionsprodlimandneighborhoodfrom the prodlim package.
