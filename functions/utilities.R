#custom asinh transform for ggplot
asinh_trans <- trans_new(
  name = "asinh",
  transform = asinh,
  inverse = sinh
)

update_geom_defaults("point", list(shape = 21,  size=1, stroke=0.3))