rndr2 <- function(x, name, ...) {
if (length(x) == 0) {
y <- clinical[[name]]
s <- rep("", length(render.default(x=y, name=name, ...)))
if (is.numeric(y)) {
p <- t.test(y ~ clinical$group)$p.value
} else {
p <- chisq.test(table(y, droplevels(clinical$group)),correct=T,simulate.p.value=T, B = 2000)$p.value
}
s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
s
} else {
render.default(x=x, name=name, ...)
}
}
