data <- read.csv("presentations/2023/RSCRUGS/bethel_algorithm.csv")


png(filename = "presentations/2023/RSCRUGS/plot0.png",
    width = 5, height = 6, units = "in", res = 500)
par(mar = c(3, 6, 0, 0))
plot(1, type = "n", xlim = c(0.5, 5.5), ylim = c(0, 0.3),
     axes = F, ann = F)
axis(side = 2, las = 1, cex = 1.5)
axis(side = 1, at = 1:5, cex = 1.5,
     labels = c("ATF", "PCod", "pollock", "POP", "SST"))
mtext(side = 2, "CV", las = 1, 4, font = 2, cex = 1.5)

points(data$CV[data$sampling.type == "SRS"],
       pch = 16, cex = 2)
points(data$CV[data$sampling.type == "single-species STRS"],
       pch = 16, col = "goldenrod", cex = 2)
legend("topleft", pch = 16, bty = "n",
       legend = c("SRS (520 stations)",
                  "SS STRS (520 stations)"),
       col = c("black", "goldenrod"))
dev.off()

png(filename = "presentations/2023/RSCRUGS/plot1.png",
    width = 5, height = 6, units = "in", res = 500)
par(mar = c(3, 6, 0, 0))
plot(1, type = "n", xlim = c(0.5, 5.5), ylim = c(0, 0.3),
     axes = F, ann = F)
axis(side = 2, las = 1, cex = 1.5)
axis(side = 1, at = 1:5, cex = 1.5,
     labels = c("ATF", "PCod", "pollock", "POP", "SST"))
mtext(side = 2, "CV", las = 1, 4, font = 2, cex = 1.5)
points(data$CV[data$sampling.type == "multi-species STRS" &
                 data$n == 422], cex = 2,
       pch = 16, col = "blue")
points(data$CV[data$sampling.type == "single-species STRS"],
       pch = 16, col = "goldenrod", cex = 2)
legend("topleft", pch = 16, bty = "n",
       legend = c("MS STRS (422 stations)",
                  "SS STRS (520 stations)"),
       col = c("blue", "goldenrod"))
dev.off()



png(filename = "presentations/2023/RSCRUGS/plot2.png",
    width = 5, height = 6, units = "in", res = 500)
par(mar = c(3, 6, 0, 0))
plot(1, type = "n", xlim = c(0.5, 5.5), ylim = c(0, 0.3),
     axes = F, ann = F)
axis(side = 2, las = 1, cex = 1.5)
axis(side = 1, at = 1:5, cex = 1.5,
     labels = c("ATF", "PCod", "pollock", "POP", "SST"))
mtext(side = 2, "CV", las = 1, 4, font = 2, cex = 1.5)
points(data$CV[data$sampling.type == "multi-species STRS" &
                 data$n == 422], cex = 2,
       pch = 16, col = "blue")
points(data$CV[data$sampling.type == "single-species STRS"],
       pch = 16, col = "goldenrod", cex = 2)

points(data$CV[data$sampling.type == "multi-species STRS" &
                 data$n == 479], cex = 2,
       pch = 16, col = "darkgrey")

legend("topleft", pch = 16, bty = "n",
       legend = c("MS STRS (422 stations)",
                  "MS STRS (479 stations)",
                  "SS STRS (520 stations)"),
       col = c("blue", "darkgrey", "goldenrod"))
dev.off()


png(filename = "presentations/2023/RSCRUGS/plot3.png",
    width = 5, height = 6, units = "in", res = 500)
par(mar = c(3, 6, 0, 0))
plot(1, type = "n", xlim = c(0.5, 5.5), ylim = c(0, 0.3),
     axes = F, ann = F)
axis(side = 2, las = 1, cex = 1.5)
axis(side = 1, at = 1:5, cex = 1.5,
     labels = c("ATF", "PCod", "pollock", "POP", "SST"))
mtext(side = 2, "CV", las = 1, 4, font = 2, cex = 1.5)
points(data$CV[data$sampling.type == "multi-species STRS" &
                 data$n == 422], cex = 2,
       pch = 16, col = "blue")
points(data$CV[data$sampling.type == "single-species STRS"],
       pch = 16, col = "goldenrod", cex = 2)

points(data$CV[data$sampling.type == "multi-species STRS" &
                 data$n == 479], cex = 2,
       pch = 16, col = "darkgrey")
points(data$CV[data$sampling.type == "multi-species STRS" &
                 data$n == 520], cex = 2,
       pch = 16, col = "black")

legend("topleft", pch = 16, bty = "n",
       legend = c("MS STRS (422 stations)",
                  "MS STRS (479 stations)",
                  "MS STRS (520 Stations)",
                  "SS STRS (520 stations)"),
       col = c("blue", "darkgrey", "black","goldenrod"))
dev.off()
