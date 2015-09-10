king.plot <- function(kingdata, gen) {

    kingdata$ped.rel <- NULL
    kingdata$ped.rel <- rep("other", length(kingdata$id.y))
    for (i in seq_along(kingdata$id.y)) {
        if (kingdata$sire.x[i] == kingdata$sire.y[i] &
            kingdata$dam.x[i] == kingdata$dam.y[i]) {
            kingdata$ped.rel[i] <- "full_sibs"
        }
        if (kingdata$id.x[i] == kingdata$sire.y[i] |
            kingdata$id.x[i] == kingdata$dam.y[i]) {
            kingdata$ped.rel[i] <- "parent_child"
        }
        if (kingdata$id.y[i] == kingdata$sire.x[i] |
            kingdata$id.y[i] == kingdata$dam.x[i]) {
            kingdata$ped.rel[i] <- "parent_child"
        }
        if (kingdata$Phi[i] < 0.05) {
            kingdata$ped.rel[i] <- "unrelated"
        }
        if (kingdata$Kinship[i] < -2) {
            kingdata$Kinship[i] <- NA
            kingdata$IBS0[i] <- NA
            kingdata$Z0[i] <- NA
        }
    }

    parent <- which(kingdata$ped.rel == "parent_child")
    sibs <- which(kingdata$ped.rel == "full_sibs")
    unrelated <- which(kingdata$ped.rel == "unrelated")
    other <- which(kingdata$ped.rel == "other")
    u <- par("usr")
    parcol <- rgb(255,0,0,80, maxColorValue = 255)
    sibcol <- rgb(255,165,0,80, maxColorValue = 255)
    othercol <- rgb(30,144,255,50,maxColorValue=255)
    unrelcol <- rgb(190,190,190,35,maxColorValue=255)

    pdf(file=paste0("IBS0_K.",gen, ".pdf"), height=5, width = 5)
    plot( kingdata$IBS0, kingdata$Kinship,
          type='p', col=rgb(190,190,190,0,maxColorValue=255), pch=".",
          main=paste0("Within-family relationship in ", gen),
          xlab="Proportion of Zero IBS",
          ylab="Genetic kinship (KING)")
    points(kingdata$IBS0[unrelated], kingdata$Kinship[unrelated],
           col=unrelcol, pch=".")
    points(kingdata$IBS0[other], kingdata$Kinship[other],
           col=othercol, pch=1)
    points(kingdata$IBS0[parent], kingdata$Kinship[parent], col=parcol, pch=4)
    points(kingdata$IBS0[sibs], kingdata$Kinship[sibs], col=sibcol, pch=6)

    abline(h = 0.3536, col = "black", lty = 3)
    abline(h = 0.1768, col = "black", lty = 3)
    abline(h = 0.0884, col = "black", lty = 3)
    abline(h = 0.0442, col = "black", lty = 3)
    abline(h = 0.0221, col = "black", lty = 3)
    legend("bottomleft", xjust=1, yjust=1, pch = c(4,6,1,20), bty="n",
           legend=c("Parent-offspring", "Full sibs", "Other", "Unrelated"),
           col = c("red", "orange", "dodgerblue", "grey"), cex=0.7)
    dev.off()



    pdf(file=paste0("PrIBD0_K.",gen,".pdf"), height=5, width = 5)
    plot( kingdata$Z0, kingdata$Kinship,
          type='p', xaxt='n', col=rgb(190,190,190,0,maxColorValue=255), pch=".",
          main=paste0("Within-family relationship in ", gen),
          xlab="Probability of Zero IBD",
          ylab="Genetic kinship (KING)")
    points(kingdata$Z0[unrelated], kingdata$Kinship[unrelated],
           col=unrelcol, pch=".")
    points(kingdata$Z0[other], kingdata$Kinship[other],
           col=othercol, pch=1)
    points(kingdata$Z0[parent], kingdata$Kinship[parent], col=parcol, pch=4)
    points(kingdata$Z0[sibs], kingdata$Kinship[sibs], col=sibcol, pch=6)

    axis(1, labels=c("0.912", "0.823", "0.646", "0.365", "0.1"), line = -0.5, lty="blank",
         at = c(.925, 0.8, 0.6464, 0.365, 0.1))
    abline(v = .1, col="black", lty=3)
    abline(v = .365, col="black", lty=3)
    abline(v = .6464, col = "black", lty=3)
    abline(v = .8232, col = "black", lty=3)
    abline(v = .9116, col = "black", lty=3)
    abline(h = 0.3536, col = "black", lty = 3)
    abline(h = 0.1768, col = "black", lty = 3)
    abline(h = 0.0884, col = "black", lty = 3)
    abline(h = 0.0442, col = "black", lty = 3)
    abline(h = 0.0221, col = "black", lty = 3)
#   u <- par("usr")
    legend("bottomleft", xjust=1, yjust=1, pch = c(4,6,1,20), bty="n",
           legend=c("Parent-offspring", "Full sibs", "Other", "Unrelated"),
           col = c("red", "orange", "dodgerblue", "grey"), cex=0.7)
    dev.off()



    pdf(file=paste0("IBS0_truncK.",gen,".pdf"), height=5, width = 5)
    plot( kingdata$IBS0, kingdata$trunK,
          type='p', col=rgb(190,190,190,0,maxColorValue=255), pch=".",
          main=paste0("Within-family relationship in ", gen),
          xlab="Proportion of Zero IBS",
          ylab="Genetic kinship (KING))")
    points(kingdata$IBS0[unrelated], kingdata$trunK[unrelated],
           col=unrelcol, pch=".")
    points(kingdata$IBS0[other], kingdata$trunK[other],
           col=othercol, pch=1)
    points(kingdata$IBS0[parent], kingdata$trunK[parent], col=parcol, pch=4)
    points(kingdata$IBS0[sibs], kingdata$trunK[sibs], col=sibcol, pch=6)

    abline(h = 0.3536, col = "black", lty = 3)
    abline(h = 0.1768, col = "black", lty = 3)
    abline(h = 0.0884, col = "black", lty = 3)
    abline(h = 0.0442, col = "black", lty = 3)
    abline(h = 0.0221, col = "black", lty = 3)
    u <- par("usr")
    legend(u[2], u[4], xjust=1, yjust=1, pch = c(4,6,1,20), bty="n",
           legend=c("Parent-offspring", "Full sibs", "Other", "Unrelated"),
           col = c("red", "orange", "dodgerblue", "grey"), cex=0.7)
    dev.off()

}

king.plot(gen50, "F50"); king.plot(gen51, "F51"); king.plot(gen52, "F52");
king.plot(gen53, "F53"); king.plot(gen54, "F54"); king.plot(gen55, "F55");
king.plot(gen56, "F56")

king.plot(gen39.40, "F39_40"); king.plot(gen40.41, "F40_41"); king.plot(gen41.42, "F41_42")
king.plot(gen42.43, "F42_43")





