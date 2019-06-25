CALAUC_train=function(time){
   
    newdata = Train[, colnames(Train) %in% c(cox_vip_list, "status", "stage","OS", "id")]
    riskScore<-predict(stepcox, type = "risk", newdata = newdata)
    risk_df <- cbind.data.frame(newdata[, c("id","OS","status")], riskScore = riskScore)
    roc <- survivalROC(Stime = risk_df$OS/365, status = risk_df$status, marker = risk_df$riskScore,predict.time = time, method = "KM")
    AUC=roc$AUC 
 }


CALAUC_test=function(time){
    newdata = Test[, colnames(Test) %in% c(cox_vip_list, "status", "stage","OS", "id")]
    riskScore<-predict(stepcox, type = "risk", newdata = newdata)
    risk_df <- cbind.data.frame(newdata[, c("id","OS","status")], riskScore = riskScore)
    roc <- survivalROC(Stime = risk_df$OS/365, status = risk_df$status, marker = risk_df$riskScore,predict.time = time, method = "KM")
    AUC=roc$AUC 
 }


 
