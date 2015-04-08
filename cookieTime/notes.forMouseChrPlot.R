### useful for making abe's "chromosome plots" without the chromosome
plot(1:100, rep(15,100), pch='|', col=rainbow(100), ylim=c(10,20))
plot(1:100, rep(15,100), pch='|', col=rainbow(100), ylim=c(10,20), cex=2)
plot(1:100, rep(15,100), pch='|', col=rainbow(100), ylim=c(10,20), cex=5)
points(1:100, rep(13,100), pch='|', col=rev(rainbow(100)), cex=5)