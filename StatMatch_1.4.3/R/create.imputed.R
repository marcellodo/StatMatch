create.imputed <- function (data.rec, data.don, mtc.ids) 
{
    # identifiers
    r.id <- mtc.ids[, 1]
    d.id <- mtc.ids[, 2]
    # records
    A <- data.rec[r.id, ]
    B <- data.don[d.id, ]
    C <- A # copy of rec
    #
    # impute NAs in a record in A
    # with values observed on the closest donor in B
    # with 2 or more NAs in a record corresponds to joint imputation
    m <- nrow(A)
    for(i in 1:m){
        tst <- is.na(A[i, ])
        C[i, tst] <- B[i, tst]
    }
    # recipient with missing imputed
    C
}