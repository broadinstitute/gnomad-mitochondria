# compile figure panel
library(ggpubr)

dir.create("figures")

source("generate_fig5A.R", echo = TRUE) 
load("figures/Fig5a.rdata")

source("generate_fig5B.R", echo = TRUE) 
load("figures/Fig5b.rdata")

source("generate_fig5C.R", echo = TRUE) 
load("figures/Fig5c.rdata")

source("generate_fig5DEF.R", echo = TRUE)
load("figures/Fig5def.rdata")

ggarrange(ggarrange(plota, plotb, labels = c("A", "B"), ncol = 2, nrow = 1, font.label = list(size = 22), widths = c(1, 1.5)),
          ggarrange(plotc, plotd, labels = c("C", "D"), ncol = 2, nrow = 1, font.label = list(size = 22), widths = c(1, 2.5)),
          ggarrange(plote, plotf, labels = c("E", "F"), ncol = 2, nrow = 1, font.label = list(size = 22), widths = c(8, 1.5)),
          ncol = 1, nrow = 3, heights = c(1, 1, 1))

#options(bitmapType = 'cairo', device = 'png')
ggsave("figures/Figure5.png", width = 25, height = 20, dpi = 300)
