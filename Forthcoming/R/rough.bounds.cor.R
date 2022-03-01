rough.bounds.cor <- function(data.A, data.B, x.lab, y.lab, z.lab, method="lm"){
    # formulas
    if(is.null(x.lab)) aa <- "."
    else aa <- paste(x.lab, collapse = "+")
    
    # 
    fyx <- paste(y.lab, aa, sep="~")
    fzx <- paste(y.lab, aa, sep="~")
    
    #### linear model
    if(method=="lm"){
        fit.y <- lm(formula = as.formula(fyx), data = data.A)
        pred.yA <- predict(fit.y)
        pred.yB <- predict(fit.y, newdata=data.B)
        
        fit.z <- lm(formula = as.formula(fzx), data = data.B)
        pred.zA <- predict(fit.z, newdata=data.A)
        pred.zB <- predict(fit.z)
    }
    #### random Forest
    if(method=="rf" | tolower(method)=="randomforest"){
        fit.y <- randomForest::randomForest(formula = as.formula(fyx), data = data.A)
        pred.yA <- predict(fit.y)
        pred.yB <- predict(fit.y, newdata=data.B)
        
        fit.z <- randomForest::randomForest(formula = as.formula(fzx), data = data.B)
        pred.zA <- predict(fit.z, newdata=data.A)
        pred.zB <- predict(fit.z)
    }
    ####### xgBoost
    if(method=="xgb" | tolower(method)=="xgboost"){
        fit.y <- caret::train(form=as.formula(fyx), data=data.A, 
                              method = "xgbTree", trControl = trainControl("cv", number = 5))
        pred.yA <- predict(fit.y)
        pred.yB <- predict(fit.y, newdata=data.B)
        
        fit.z <- caret::train(form=as.formula(fzx), data=data.B, 
                              method = "xgbTree", trControl = trainControl("cv", number = 5))
        pred.zA <- predict(fit.z, newdata=data.A)
        pred.zB <- predict(fit.z)
        
    }
    df.A <- data.frame(data.A[,y.lab], pred.y=pred.yA, pred.z=pred.zA)
    df.B <- data.frame(pred.y=pred.yB, pred.z=pred.zB, data.B[,z.lab])
    
        
}