GGsurvplot = function(lncrna_list,dataset,stepcox) { 
      
	newdata = dataset[, colnames(dataset) %in% c(lncrna_list, "status", "stage","OS", "id")]
      riskScore <-predict(stepcox, type = "risk", newdata = newdata)
      risk <- as.factor(ifelse(riskScore > median(riskScore), "high", "low"))
      risk_df <- cbind.data.frame(newdata[, c("id","OS","status")], riskScore = riskScore, risk = risk)
      risk_df <- risk_df%>% mutate(status = as.numeric(as.vector(status)))
      fit <- survfit(Surv(OS/365, status) ~ risk, data = risk_df)
      g=ggsurvplot(fit, data = risk_df, risk.table = F, pval = T, ggtheme = theme_survminer(),
	  title = "Overall survival probability of TCGA-COAD in predicted risk groups",
	  xlab = "Time(years)", legend = "top", surv.median.line = "hv")
      s=survdiff(formula = Surv(OS/365, status) ~ risk, data = risk_df)
	  return(c(g,s))
}
