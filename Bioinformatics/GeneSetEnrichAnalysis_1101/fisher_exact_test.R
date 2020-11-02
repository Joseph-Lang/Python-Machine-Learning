library('MASS')
head(Cars93)

tab1 <- table(Cars93$Type, Cars93$Origin)
# Chi-squared test
chisq.test(Cars93$Type, Cars93$Origin)
# Fisher's exact test
fisher.test(Cars93$Type, Cars93$Origin)

p.val <- matrix(c(3, 297, 37, 29663), nrow = 2)
rownames(p.val) <- c("In Pathway", "Not In Pathway")
colnames(p.val) <- c("User Genes", "Genome")
# Fisher's exact test
fisher.test(p.val)
# EASE score
score <- matrix(c(2, 297, 38, 29663), nrow = 2)
rownames(score) <- c("In Pathway", "Not In Pathway")
colnames(score) <- c("User Genes", "Genome")
fisher.test(score)
