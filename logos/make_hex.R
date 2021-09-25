library(hexSticker)
# Source image
img <- "logos/jSDM_wide.png"
# Greyscale version
sticker(subplot=img, s_width=0.58, s_x=1, s_y=0.72,
        package="jSDM", p_size=30, p_x=1, p_y=1.4,
        h_fill=grey(0.7), h_color=grey(0.4),
        h_size=1.3,
        url="https://ecology.ghislainv.fr/jSDM",
        u_color="black",
        u_size=3,
        filename="logos/logo_grey.png")
# Colored version
# https://www.color-hex.com/color-palette/78245
sticker(subplot=img, s_width=0.58, s_x=1, s_y=0.72,
        package="jSDM", p_size=30, p_x=1, p_y=1.4,
        p_color="#bd4545",
        h_fill="#f0d171", h_color="#9a542e",
        h_size=1.3,
        url="https://ecology.ghislainv.fr/jSDM",
        u_color="#6d1950",
        u_size=3,
        filename="man/figures/logo.png")