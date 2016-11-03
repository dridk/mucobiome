data = read.table("all_report.txt", as.is=T , header=T)
barplot(t(as.matrix(data[order(data$raw),c("raw","merge","clean")])), 
        names=data$sample,
        las=2,
        cex.names=0.5,
        width=0.1, 
        main="Reads filtering per sample",
        col = grey.colors(3)
        )

legend("topleft", c("raw","merge","clean"), fill =grey.colors(3))