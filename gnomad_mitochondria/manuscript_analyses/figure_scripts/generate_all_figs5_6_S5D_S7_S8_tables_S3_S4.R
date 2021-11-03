# do analyses and make figures related to patterns of variation, and pathogenic variants
dir.create("figures")
dir.create("tables")

# main text figures
source("collate_figure5.R", echo = TRUE) 
source("generate_fig6.R", echo = TRUE) 

# supplementary figures
source("generate_figS5D.R", echo = TRUE) 
source("generate_figS7.R", echo = TRUE) 
source("generate_figS8.R", echo = TRUE)
ggsave("figures/FigureS8.png", width = 15, height = 13) #due to a viewport error it isn't saving the figure within the script, so providing command here

# supplementary tables
source("generate_tableS3.R", echo = TRUE) 
source("generate_tableS4.R", echo = TRUE) 





