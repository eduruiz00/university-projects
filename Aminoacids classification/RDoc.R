#----
#librerias
#install.packages(c("ChemometricsWithR", "MASS"))
#install.packages(c("e1071", "sfsmisc", "class", "caret"))
#ChemometricsWithR::installChemometricsWithRData()
library("stats")
library("ChemometricsWithR")
library("MASS")
library("pls")
library("sfsmisc")
library("e1071")
library("class")
library("lattice")
library("ggplot2")
library("caret")
set.seed(27)
#----
#load data
x1 <- as.matrix(read.table("cancer.dat"))
x2 <- as.matrix(read.table("smoker.dat"))
x3 <- as.matrix(read.table("healthy.dat"))

#----
#classification
num_samples_per_class_in_train <- round(0.75*200)
num_idx <- 1:200

cancer_idx_train <- sample(num_idx, num_samples_per_class_in_train)
cancer_idx_test <- setdiff(num_idx, cancer_idx_train)

smoker_idx_train <- sample(num_idx, num_samples_per_class_in_train)
smoker_idx_test <- setdiff(num_idx, smoker_idx_train)

health_idx_train <- sample(num_idx, num_samples_per_class_in_train)
health_idx_test <- setdiff(num_idx, health_idx_train)

# use the indexes to split the matrix into train and test
x1_t <- x1[cancer_idx_train,]
x2_t <- x2[smoker_idx_train,]
x3_t <- x3[health_idx_train,]
x1_v <- x1[cancer_idx_test,]
x2_v <- x2[smoker_idx_test,]
x3_v <- x3[health_idx_test,]

xt <- rbind(x1_t,x2_t,x3_t)
xv <- rbind(x1_v,x2_v,x3_v)
clt <- c(rep(1,150),rep(2,150),rep(3,150))
clv <- c(rep(1,50),rep(2,50),rep(3,50))
total <- rbind(xt,xv)

#----
#hacemos el scatter plot
colours_legend <- c("1"=2, "2"=3, "3"=4)
point_shapes_legend <- c(train = 0, validation = 8)
point_shapes <- point_shapes_legend[c(rep("train", nrow(xt)),
                                      rep("validation", nrow(xv)))]
colours <- colours_legend[c(clt, clv)]
pairs(total[,1:5], labels = paste("AA", 1:5), pch = point_shapes,
      col = colours)

#----
#k-nearest
k=11
subject_condition_tst_pca_knn_pred <- class::knn(train = xt, test = xv,
                                                 cl = clt, k = k)

confmat_pca_knn <- table(clv, subject_condition_tst_pca_knn_pred)
print(confmat_pca_knn)

CR_pca_knn <- sum(diag(confmat_pca_knn))/sum(confmat_pca_knn)
message("The classification rate for kNN is: ", 100*round(CR_pca_knn, 2), "%")
#----
#linear discriminant classifier
lda_model <- lda(xt, grouping = clt)
lda_pred_test <- predict(lda_model, xv)
lda_scores_tst <- lda_pred_test$x
subject_condition_tst_lda_pred <- lda_pred_test$class
lda_pred_trn <- predict(lda_model, xt)
lda_scores_trn <- lda_pred_trn$x

plot(x = lda_scores_trn[, 1],
     y = lda_scores_trn[, 2],
     col = colours_legend[clt],
     pch = point_shapes_legend["train"],
     xlab = "LD1", ylab = "LD2", main = "LDA scores before PCA",
     cex.main=1, cex.lab=0.8, cex.axis=0.8)
points(x = lda_scores_tst[, 1],
       y = lda_scores_tst[, 2],
       col = colours_legend[clv],
       pch = point_shapes_legend["train"])
legend("bottomright",
       legend = c("Lung cancer","Smoker", "Healthy"),
       pch = 19,
       col = colours_legend,
       title = "Condition",
       cex=0.8)
legend("bottomleft",
       legend = c("Train","Validation"), 
       pch = point_shapes_legend,
       title = "Subset")

confmat_lda <- table(clv, subject_condition_tst_lda_pred)
print(confmat_lda)

CR_lda <- sum(diag(confmat_lda))/sum(confmat_lda)
message("The classification rate for LDA is ", 100*round(CR_lda, 2), "%")

#----
#PCA reduction
#mean_spectrum_trn <- colMeans(xt)
#intensity_trn_preproc <- matrix(0, nrow = nrow(xt), ncol = ncol(xt))
#for (i in 1:ncol(intensity_trn_preproc)) {
#   intensity_trn_preproc[, i] <- xt[, i] - mean_spectrum_trn[i] }

#intensity_tst_preproc <- matrix(0, nrow = nrow(xv), ncol = ncol(xv))
#intensity_tst_preproc[, i] <- xv[, i] - mean_spectrum_trn[i] }

intensity_trn_preproc <- scale(xt, center=TRUE, scale=FALSE)
computed_mean <- attr(intensity_trn_preproc, "scaled:center") 
intensity_tst_preproc <- scale(xv, center = computed_mean, scale = FALSE)

pca_model <- PCA(intensity_trn_preproc)
summary(pca_model)

pca_model_var <- variances(pca_model)
pca_model_var_percent <- 100*cumsum(pca_model_var)/sum(pca_model_var)
plot(x = 1:length(pca_model_var_percent),
     y = pca_model_var_percent,
     type = "l",
     xlab = "Number of principal components", ylab = "Cummulated variance (%)",
     main = "% variance captured vs #PC used")
pca_model_loadings <- loadings(pca_model, 3)
matplot(pca_model_loadings,col=c(1,4,3),type="l", xlab="Aminoacids",
        main="Loadings for the 3 principal components",ylab="Loadings",lty=1)
legend("bottomright",
       legend = c("PC1", "PC2", "PC3"),
       pch = 19,
       title = "Condition",
       col = c(1,4,3),
       cex=0.6)
abline(v=c(2,12,20),col=2)
#----
#SCATTER PLOT PCA
pca_model_ncomp <- 3
intensity_scores_trn <- scores(pca_model, npc = pca_model_ncomp)
intensity_scores_test <- project(pca_model, npc = pca_model_ncomp,
                                 intensity_tst_preproc)


intensity_scores_all <- rbind(intensity_scores_trn, intensity_scores_test)
colours <- colours_legend[c(clt, clv)]
pairs(intensity_scores_all[,1:3], labels = paste("AA", 1:3), pch = point_shapes,
      col = colours)
#----
plot(x = intensity_scores_all[,1], y = intensity_scores_all[,2], pch = point_shapes,
     col = colours,
     main = "PCA Score plot",
     xlab = "PC 1",
     ylab = "PC 2")
legend("bottomright",
       legend = c("Lung cancer", "Smoker", "Healthy"),
       pch = 19,
       title = "Condition",
       col = colours_legend)
legend("topright",
       legend = c("Train","Validation"),
       pch = point_shapes_legend,
       title = "Subset",
       col = "black")

#----
#pca-knn
subject_condition_tst_pca_knn_pred <- class::knn(train = intensity_scores_trn,
                                                 test = intensity_scores_test,
                                                 cl = clt, k = 11)
confmat_pca_knn <- table(clv, subject_condition_tst_pca_knn_pred)
print(confmat_pca_knn)

CR_pca_knn <- sum(diag(confmat_pca_knn))/sum(confmat_pca_knn)
message("The classification rate for PCA-kNN is: ", 100*round(CR_pca_knn, 2), "%")

#----
#PCA-LDA
pca_lda_model <- lda(intensity_scores_trn, grouping = clt)
pca_lda_pred_test <- predict(pca_lda_model, intensity_scores_test)

lda_scores_tst1 <- pca_lda_pred_test$x
subject_condition_tst_pca_lda_pred <- pca_lda_pred_test$class

pca_lda_pred_trn <- predict(pca_lda_model, intensity_scores_trn)
lda_scores_trn1 <- pca_lda_pred_trn$x
plot(x = lda_scores_trn1[, 1],
     y = lda_scores_trn1[, 2],
     col = colours_legend[clt],
     pch = point_shapes_legend["train"],
     xlab = "LD1", ylab = "LD2", main = "LDA scores, after PCA reduction")
points(x = lda_scores_tst1[, 1],
       y = lda_scores_tst1[, 2],
       col = colours_legend[clv],
       pch = point_shapes_legend["validation"])
legend("bottomright",
       legend = c("Lung cancer","Smoker","Healthy"),
       pch = 19,
       col = colours_legend,
       title = "Condition")
legend("bottomleft",
       legend = c("Train","Validation"),
       pch = point_shapes_legend,
       title = "Subset")

confmat_pca_lda <- table(clv, subject_condition_tst_pca_lda_pred)
print(confmat_pca_lda)

CR_pca_lda <- sum(diag(confmat_pca_lda))/sum(confmat_pca_lda)
message("The classification rate for PCA-LDA is ", 100*round(CR_pca_lda, 2), "%")

