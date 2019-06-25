### Give a cox_vip_list years, Train  or Validate data set(type: data frame)and a coxph model for perdiction ,the default year is 5.


PlotsurvROC = function(predict.time=5,Train,Validate, model = stepcox) {
	   require(survminer)
	   require(survivalROC)
      
      newdata_t = Train[, colnames(Train) %in% c(cox_vip_list, "status", "stage","OS", "id")]
      riskScore_t <-predict(stepcox, type = "risk", newdata = newdata_t)
      risk_t <- as.factor(ifelse(riskScore_t > median(riskScore_t), "high", "low"))
      risk_df_t <- cbind.data.frame(newdata_t[, c("id","OS","status")], riskScore = riskScore_t, risk = risk_t)
      
      par(omi = c(0.5, 0.5, 0.5, 0.5), font.lab = 1, font.axis = 1)
      roc_t <- survivalROC(Stime = risk_df_t$OS/365, status = risk_df_t$status, marker = risk_df_t$riskScore,
				   predict.time = predict.time, method = "KM")##
	  plot(roc_t$FP, roc_t$TP, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "#FC4E07",
      xlab = "False positive rate", ylab = "True positive rate",
	  main = paste("ROC of multivariable CoxPH model for risk prediction\n (Method = KM,Year =", predict.time,")"),
	  lwd = 2, cex.main = 1, cex.lab = 1, cex.axis = 1, font = 1, las = 1)
	  

     with(roc_t, points(FP[with(roc_t, which.min(1 - TP + FP))], TP[with(roc_t, which.min(1 - TP + FP))], cex = 2, col = "#FC4E07"))
     with(roc_t, text(FP[with(roc_t, which.min(1 - TP + FP))], TP[with(roc_t, which.min(1 - TP + FP))], cex = 1, col = "black", labels = paste(round(with(roc_t, min(1 - TP + FP)),3)), xpd = T))
     ################
      newdata_v = Validate[, colnames(Validate) %in% c(cox_vip_list, "status", "stage","OS", "id")]
      riskScore_v <-predict(stepcox, type = "risk", newdata = newdata_v)
      risk_v <- as.factor(ifelse(riskScore_v > median(riskScore_v), "high", "low"))
      risk_df_v <- cbind.data.frame(newdata_v[, c("id","OS","status")], riskScore = riskScore_v, risk = risk_v)
      
      roc_v <- survivalROC(Stime = risk_df_v$OS/365, status = risk_df_v$status, marker = risk_df_v$riskScore,
				   predict.time = predict.time, method = "KM")##
	  
      lines(roc_v$FP, roc_v$TP, type = "l", col = "#0073C2FF",lwd = 2)
	  with(roc_v, which.min(1 - TP + FP))
      with(roc_v, min(1 - TP + FP))
      #
      with(roc_v, points(FP[with(roc_v, which.min(1 - TP + FP))], TP[with(roc_v, which.min(1 - TP + FP))], cex = 2, col = "#0073C2FF"))
      with(roc_v, text(FP[with(roc_v, which.min(1 - TP + FP))], TP[with(roc_v, which.min(1 - TP + FP))], cex = 1, col = "black", labels = paste(round(with(roc_v, min(1 - TP + FP)),3)), xpd = T))
	  abline(0, 1)
	 
	 legend("bottomright", col = c("#FC4E07", "#0073C2FF"), bty = "n",
	 cex = 0.8, lwd = 2, title = "Model performace ", title.adj = 0.5,
	 legend = c(paste("Training data AUC", "(", round(roc_t$AUC, 3), ")"),
	 paste("Testing data AUC", "(", round(roc_v$AUC, 3), ")")))
	}

