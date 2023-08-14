library(dartR)

res_fstp <- pairwise_table(x=platypus.gl,
                      pw_fun = FSTp)


res_detp <- pairwise_table(x=platypus.gl,
                                   pw_fun = Dest)
