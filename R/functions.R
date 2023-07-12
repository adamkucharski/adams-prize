## little man function

draw_man <- function(x0 = 0, y0 = 0, radius = 0.2, height = 1, col = "grey", lwd = 2)
{
  draw.circle(x = x0+radius, y = y0+height-radius, radius = radius, border = col, lwd = lwd)
  segments(x0+radius, y0+radius, x0+radius, y0+height-2*radius, col = col, lwd = lwd)
  segments(x0, y0, x0+radius, y0+radius, col = col, lwd = lwd)
  segments(x0+2*radius, y0, x0+radius, y0+radius, col = col, lwd = lwd)
  segments(x0, y0+1.5*radius, x0+radius, y0+2.5*radius, col = col, lwd = lwd)
  segments(x0+2*radius, y0+1.5*radius, x0+radius, y0+2.5*radius, col = col, lwd = lwd)
}
