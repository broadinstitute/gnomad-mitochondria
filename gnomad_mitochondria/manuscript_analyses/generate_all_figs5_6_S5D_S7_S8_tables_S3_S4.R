# do analyses and make figures related to patterns of variation, and pathogenic variants
dir.create("figures")
dir.create("tables")

# main text figures
source("collate_Figure5.R", echo = TRUE) 
source("generate_Fig6.R", echo = TRUE) 

# supplementary figures
source("generate_FigS5D.R", echo = TRUE) 
source("generate_FigS7.R", echo = TRUE) 
source("generate_FigS8.R", echo = TRUE)
ggsave("figures/FigureS8.png", width = 15, height = 13) #due to a viewport error it isn't saving the figure within the script, so providing command here

# supplementary tables
source("generate_tableS3.R", echo = TRUE) 
source("generate_tableS4.R", echo = TRUE) 





