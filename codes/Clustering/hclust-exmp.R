set.seed(0)
x = rbind(scale(matrix(rnorm(2*20),ncol=2),cent=c(1,1),scale=F),
          scale(matrix(rnorm(2*30),ncol=2),cent=-c(1,1),scale=F))
x = rbind(x,matrix(runif(2*10,min(x),max(x)),ncol=2))
d = dist(x)

plot(x)

tree.sing = hclust(d,method="single")
tree.comp = hclust(d,method="complete")
tree.avg = hclust(d,method="average")

par(mfrow=c(2,3))
plot(tree.sing,labels=F,hang=-1e-10)
plot(tree.comp,labels=F,hang=-1e-10)
plot(tree.avg,labels=F,hang=-1e-10)

labs.sing = cutree(tree.sing,k=3)
labs.comp = cutree(tree.comp,k=3)
labs.avg = cutree(tree.avg,k=3)

cols = c("red","darkgreen","blue")
plot(x,col=cols[labs.sing])
plot(x,col=cols[labs.comp])
plot(x,col=cols[labs.avg])

# Compare clustering with 3 and 4 groups
cols = c("red","blue","darkgreen","purple")
labs.avg2 = cutree(tree.avg,k=4)

par(mfrow=c(1,2))
plot(x,col=cols[labs.avg],main="3 clusters")
plot(x,col=cols[labs.avg2],main="4 clusters")

table(labs.avg,labs.avg2)