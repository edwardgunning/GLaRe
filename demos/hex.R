install.packages("hexSticker")
library(hexSticker)
sticker("inst/figures/sticker-image.png",
        package="GlaRE",
        p_size=25,
        s_x=1,
        s_y=.85,
        s_width=.7525,
        p_color = "black",
        h_color = "black",
        h_fill = "white",
        filename="inst/figures/hex.png")
