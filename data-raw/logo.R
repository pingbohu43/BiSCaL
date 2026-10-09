# Creates man/figures/logo.png, the package logo (run from the package root):
# two chains of responses over visits with cross-lagged arrows.
png("man/figures/logo.png", width = 518, height = 600, bg = "transparent", res = 144)
op <- par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(-1, 1), ylim = c(-1.155, 1.155), asp = 1, axes = FALSE, xlab = "", ylab = "")
ang <- seq(90, 450, by = 60) * pi / 180
polygon(1.1 * cos(ang), 1.1 * sin(ang), col = "#1f4e5f", border = "#2a78d6", lwd = 7)
xs <- c(-0.55, -0.05, 0.45)
y_top <- 0.38
y_bot <- -0.08
# horizontal links (each response depends on its previous value)
arrows(xs[-3] + 0.12, y_top, xs[-1] - 0.12, y_top, col = "white", lwd = 4, length = 0.07)
arrows(xs[-3] + 0.12, y_bot, xs[-1] - 0.12, y_bot, col = "white", lwd = 4, length = 0.07)
# cross links: Y2[k-1] -> Y1[k] and Y1[k] -> Y2[k]
arrows(xs[-3] + 0.1, y_top - 0.07, xs[-1] - 0.1, y_bot + 0.07, col = "#9cc3ee", lwd = 3,
       length = 0.06)
arrows(xs[-1], y_bot + 0.12, xs[-1], y_top - 0.12, col = "#9cc3ee", lwd = 3, length = 0.06)
points(xs, rep(y_top, 3), pch = 21, bg = "#eb6834", col = "white", cex = 3.2, lwd = 3)
points(xs, rep(y_bot, 3), pch = 21, bg = "#2a78d6", col = "white", cex = 3.2, lwd = 3)
text(0, -0.62, "BiSCaL", col = "white", cex = 1.75, font = 2)
par(op)
dev.off()
