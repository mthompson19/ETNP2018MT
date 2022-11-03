# Oxygen Depth Profiles 
## Madeleine Thompson
### ETNP 16s data 

library("ggplot2"); packageVersion("ggplot2")
library(ggpubr)

##############


oxygenstn1data <- read.csv("data/Station1oxy.csv")

oxygenstn2data <- read.csv("data/station2oxy.csv")

oxygenstn3data <- read.csv("data/station3oxy.csv")

#Depth Profiles 
oxydepth1 <- ggplot2::ggplot(oxygenstn1data,aes(x=depth..m.,y=oxygen..uM.))+
  geom_line(size=1)+
  geom_point()+
  scale_x_reverse()+
  scale_y_continuous(position="right")+
  coord_flip() +
  theme_classic() +
  ylab(expression(paste("Oxygen Concentration (", mu, "M)"))) +
  xlab("Depth (m)") +
  ggtitle("Station 1") +
  theme(text=element_text(size=20))
oxydepth1

oxydepth2 <- ggplot2::ggplot(oxygenstn2data,aes(x=depth..m.,y=oxygen..uM.))+
  geom_line(size=1)+
  geom_point()+
  scale_x_reverse()+
  scale_y_continuous(position="right")+
  coord_flip() +
  theme_classic() +
  ylab(expression(paste("Oxygen Concentration (", mu, "M)"))) +
  xlab("Depth (m)") +
  ggtitle("Station 2") +
  theme(text=element_text(size=20))
oxydepth2

oxydepth3 <- ggplot2::ggplot(oxygenstn3data,aes(x=depth..m.,y=oxygen..uM.))+
  geom_line(size=1)+
  geom_point()+
  scale_x_reverse()+
  scale_y_continuous(position="right")+
  coord_flip() +
  theme_classic() +
  ylab(expression(paste("Oxygen Concentration (", mu, "M)"))) +
  xlab("Depth (m)") +
  ggtitle("Station 3") +
  theme(text=element_text(size=20))
oxydepth3

depthfigure <- ggpubr::ggarrange(oxydepth1, oxydepth2, oxydepth3, 
                                 labels = c("A", "B", "C"), 
                                 font.label = list(size = 25, color = "black", face = "bold", family = NULL),
                                 ncol = 3, nrow = 1)
ggplot2::ggsave("figures/profiles.png", plot=depthfigure, device="png",
                scale=1, width = 40, height=25, units=c("cm"), dpi=300, limitsize = FALSE)


############

